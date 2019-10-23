/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"
#include <ctime>
#include <chrono>

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, channels);
    REQ_ARG(std::string, period);
    REQ_ARG(std::string, selections);
    REQ_ARG(std::string, unc_sources);
    REQ_ARG(std::string, jet_orderings);
    REQ_ARG(bool, isData);
    REQ_ARG(bool, hasBjetPair);
    OPT_ARG(bool, runSVFit, true);
    OPT_ARG(bool, runKinFit, true);
    OPT_ARG(Long64_t, max_events_per_tree, std::numeric_limits<Long64_t>::max());
    OPT_ARG(std::string, working_path, "./");
};

namespace analysis {

class DummyCacheProducer_t {
public:
    using CacheEvent = cache_tuple::CacheEvent;
    using CacheTuple = cache_tuple::CacheTuple;
    using CacheSummaryTuple = cache_ntuple::CacheSummaryTuple;
    using clock = std::chrono::system_clock;

    DummyCacheProducer_t(const Arguments& _args) : args(_args), outputFile(root_ext::CreateRootFile(args.output_file())),
                cacheSummary("summary", outputFile.get(), false), start(clock::now()), run_period(Parse<analysis::Period>(args.period())),
                progressReporter(10, std::cout)
    {
        EventCandidate::InitializeJecUncertainty(run_period,args.working_path());

        auto signalModes = SplitValueListT<analysis::SignalMode>(args.selections(),false,",");
        for(unsigned n = 0; n < signalModes.size(); ++n){
            signalObjectSelectors.emplace_back(signalModes.at(n));
        }
        unc_sources = SplitValueListT<analysis::UncertaintySource>(args.unc_sources(),false,",");
        vector_jet_ordering = SplitValueListT<JetOrdering>(args.jet_orderings(),false,",");
        channels = SplitValueList(args.channels(),false,",");

        cacheSummary().numberOfOriginalEvents = 0;
        cacheSummary().numberOfTimesSVFit = 0;
        cacheSummary().numberOfTimesKinFit = 0;

    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% selection.\n")
                   % args.input_file() % args.output_file() % args.selections();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        size_t n_tot_events = 0;
        std::map<std::string,std::pair<std::shared_ptr<ntuple::EventTuple>,Long64_t>> map_event;
        for(unsigned c = 0; c < channels.size(); ++c){
            auto originalTuple = ntuple::CreateEventTuple(channels.at(c),originalFile.get(),true,ntuple::TreeState::Full);
            const Long64_t n_events = std::min(args.max_events_per_tree(),originalTuple->GetEntries());
            map_event[channels.at(c)] = std::make_pair(originalTuple,n_events);
            n_tot_events += n_events;
        }

        size_t n_processed_events = 0;
        progressReporter.SetTotalNumberOfEvents(n_tot_events);
        for(unsigned c = 0; c < channels.size(); ++c){
            std::cout << "Channel: " << channels.at(c) << std::endl;
            CacheTuple cache(channels.at(c), outputFile.get(), false);
            auto& originalTuple = *map_event.at(channels.at(c)).first;
            const Long64_t n_entries = map_event.at(channels.at(c)).second;
            for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
                originalTuple.GetEntry(current_entry);
                if(static_cast<Channel>(originalTuple().channelId) == Channel::MuMu){ //temporary fix due tue a bug in mumu channel in production
                    originalTuple().first_daughter_indexes = {0};
                    originalTuple().second_daughter_indexes = {1};
                }
                FillCacheTuple(cache, originalTuple.data());
                cache.Fill();
                if(++n_processed_events % 100 == 0) progressReporter.Report(n_processed_events,false);
            }
            progressReporter.Report(n_processed_events, true);
            cache.Write();
            const auto stop = clock::now();
            cacheSummary().exeTime = static_cast<UInt_t>(std::chrono::duration_cast<std::chrono::seconds>(stop - start).count());
            cacheSummary.Fill();
            cacheSummary.Write();
        }
        progressReporter.Report(n_tot_events,true);




    }

private:

    void FillCacheTuple(CacheTuple& cacheTuple, const ntuple::Event& event)
    {
        std::set<size_t> Htt_indexes;
        std::set<std::pair<size_t,size_t>> HH_indexes;
        for(unsigned selector = 0; selector < signalObjectSelectors.size(); ++selector){
            for(unsigned ordering = 0; ordering < vector_jet_ordering.size(); ++ordering){
                for(unsigned source = 0; source < unc_sources.size(); ++source){
                    for(int variation = -1; variation < 2; ++variation){
                        cacheSummary().numberOfOriginalEvents++;
                        cacheTuple().run = event.run;
                        cacheTuple().lumi = event.lumi;
                        cacheTuple().evt = event.evt;
                        const SignalObjectSelector& signalObjectSelector = signalObjectSelectors.at(selector);
                        JetOrdering jet_ordering = vector_jet_ordering.at(ordering);
                        UncertaintySource unc_source = unc_sources.at(source);
                        UncertaintyScale scale = static_cast<UncertaintyScale>(variation);
                        if(scale != UncertaintyScale::Central && unc_source == UncertaintySource::None) continue;
                        if(scale == UncertaintyScale::Central && unc_source != UncertaintySource::None) continue;

                        std::cout << "event.jets_p4.size(): " << event.jets_p4.size() << std::endl;
                        for(unsigned n = 0; n < event.jets_p4.size(); ++n){
                            cacheTuple().jet_hh_score_index.push_back(n);
                            cacheTuple().jet_hh_score_unc_scale.push_back(variation);
                            cacheTuple().jet_hh_score_unc_source.push_back(static_cast<Int_t>(unc_sources.at(source)));
                            cacheTuple().jet_hh_score_index.push_back(1.0);
                        }

                    }
                }
            }

        }

    }

private:
    Arguments args;
    std::shared_ptr<TFile> outputFile;
    CacheSummaryTuple cacheSummary;
    const clock::time_point start;
    analysis::Period run_period;
    std::vector<std::string> channels;
    std::vector<SignalObjectSelector> signalObjectSelectors;
    std::vector<UncertaintySource> unc_sources;
    std::vector<JetOrdering> vector_jet_ordering;
    analysis::tools::ProgressReporter progressReporter;
};

} // namespace analysis

PROGRAM_MAIN(analysis::DummyCacheProducer_t, Arguments)
