/*! Analyzer
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */


#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"
#include "hh-bbtautau/Analysis/source/ProcessAnaTuple.cxx"
#include <string>
#include <iostream>
#include <unordered_map>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <TSystemDirectory.h>
#include <TSystemFile.h>

namespace analysis {

struct AnalyzerArguments : CoreAnalyzerArguments {
    REQ_ARG(Channel, channel);
    REQ_ARG(std::string, input);
    REQ_ARG(std::string, signal_name);
    REQ_ARG(std::string, output);
    OPT_ARG(std::string, var, "");
};

class EvaluateTauUncertainties : public EventAnalyzerCore {
public:
    using AnaData = ::analysis::EventAnalyzerData;
    using AnaDataCollection = ::analysis::EventAnalyzerDataCollection;
    using PlotsProducer = ::analysis::StackedPlotsProducer;

    EvaluateTauUncertainties(const AnalyzerArguments& _args) :
        EventAnalyzerCore(_args, _args.channel()), args(_args), activeVariables({args.var()}),
        outputFile(root_ext::CreateRootFile(args.output()))
    {
        histConfig.Parse(FullPath(ana_setup.hist_cfg));
    }

    void Run()
    {
        const std::set<std::string> signal_names(ana_setup.signals.begin(), ana_setup.signals.end());
        const std::set<std::string> bkg_names(ana_setup.backgrounds.begin(), ana_setup.backgrounds.end());
        const std::set<std::string> data_names(ana_setup.data.begin(), ana_setup.data.end());
        std::set<std::string> all_samples;
        all_samples.insert(signal_names.begin(),signal_names.end());
        all_samples.insert(bkg_names.begin(),bkg_names.end());
        all_samples.insert(data_names.begin(),data_names.end());

        const auto samplesToDraw = PlotsProducer::CreateOrderedSampleCollection(
                    ana_setup.draw_sequence, sample_descriptors, cmb_sample_descriptors, ana_setup.signals,
                    ana_setup.data, args.channel());

        std::map<analysis::EventCategory,size_t> categories_map = {{ analysis::EventCategory::TwoJets_OneBtag_Resolved_noVBF(), 0 },
                                                                   { analysis::EventCategory::TwoJets_TwoBtagPlus_Resolved_noVBF(), 1 },
                                                                   { analysis::EventCategory::TwoJets_TwoLooseBtagPlus_Boosted_noVBF(), 2},
                                                                   { analysis::EventCategory::FourJets_OneBtagPlus_VBF(), 3}};

        auto inputFile = root_ext::OpenRootFile(args.input());
        AnaDataCollection anaDataCollection(outputFile, channelId, nullptr, activeVariables,
                                            histConfig.GetItems(), all_samples, nullptr);

        for(const auto& category : ana_setup.categories){
            std::cout<<category<<std::endl;
            if(!categories_map.count(category))
                continue;
            size_t n = categories_map.at(category);
            for(const auto& sample_name : all_samples){
                SampleDescriptor& sample = sample_descriptors.at(sample_name);
                for(const SampleDescriptorBase::Point& item : sample.working_points){
                    for (const auto& subcategory : sub_categories_to_process){

                        const EventAnalyzerDataId anaDataId(category, subcategory, analysis::EventRegion::OS_Isolated(),
                                                            analysis::EventEnergyScale::Central, item.full_name);
                        if(item.datacard_name.empty()) continue;

                        if(ana_setup.IsSignal(sample.name) && item.full_name != args.signal_name())
                            continue;

                        const std::string hist_dir_name = ana_setup.IsSignal(sample.name) ?
                                "hh_ttbb_"+ToString(args.channel())+"_"+ToString(n)+"_13TeV_postfit/"+sample.postfit_name :
                                "hh_ttbb_"+ToString(args.channel())+"_"+ToString(n)+"_13TeV_postfit/"+item.datacard_name ;

                        auto hist = std::shared_ptr<TH1>
                                (root_ext::TryReadObject<TH1>(*inputFile,hist_dir_name));
                        if (!hist){
                            std::cout << "WARNING! Datacard Name doesn't correspond to sample name: " << anaDataId <<
                                         ", hist dir name: "<< hist_dir_name <<  std::endl;
                            continue;
                        }
                        std::cout << "Loading " << anaDataId << " ..." << std::endl;
                        auto& histAnaData = anaDataCollection.Get(anaDataId);
                        auto& histAnaData_entry = histAnaData.GetEntryEx<TH1D>(args.var());
                        histAnaData_entry().CopyContent(*hist);

                    }
                }
            }
        }

        std::cout << std::endl;
        std::cout << "\tProcessing combined samples " << std::endl;

        for(const auto& subcategory : sub_categories_to_process) {
            ProcessCombinedSamples(anaDataCollection, subcategory, ana_setup.cmb_samples);
        }

        std::cout << "\t\tCreating plots..." << std::endl;
        PlotsProducer plotsProducer(anaDataCollection, samplesToDraw, FullPath(ana_setup.plot_cfg),
                                    ana_setup.plot_page_opt);
        std::string pdf_prefix = args.output();
        plotsProducer.PrintStackedPlots(pdf_prefix, EventRegion::SignalRegion(), ana_setup.categories,
                                        sub_categories_to_process, signal_names, &sample_descriptors.at("TotalBkg"));

    }

//private:



private:
    AnalyzerArguments args;
    std::set<std::string> activeVariables;
    std::shared_ptr<TFile> outputFile;
    PropertyConfigReader histConfig;
};

} // namespace analysis

PROGRAM_MAIN(analysis::EvaluateTauUncertainties, Arguments) // definition of the main program function
