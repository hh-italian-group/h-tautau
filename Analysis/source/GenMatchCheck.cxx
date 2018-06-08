/*! Analyzer to check genMatch information.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/SyncTupleHTT.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "AnalysisTools/Core/include/AnalyzerData.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, period);
};

namespace analysis {

class GenMatchCheckData : public root_ext::AnalyzerData {
public:
    using AnalyzerData::AnalyzerData;

    TH1D_ENTRY(mass_bothGood, 50, 0, 300)
    TH1D_ENTRY(mass_bothBad, 50, 0, 300)
    TH1D_ENTRY(mass_1stGood_2ndBad, 50, 0, 300)
    TH1D_ENTRY(mass_1stBad_2ndGood, 50, 0, 300)
    TH1D_ENTRY(mass_allBad, 50, 0, 300)

};

class GenMatchCheck {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;

    GenMatchCheck(const Arguments& _args) : args(_args), output(root_ext::CreateRootFile(args.output_file())),
        anaData(output) {
            run_period = analysis::EnumNameMap<analysis::Period>::GetDefault().Parse(args.period());
        }

    void Run()
    {
        static const std::map<Channel, std::pair<Float_t,Float_t>> genMatchInfo = {
            { Channel::ETau, { 3,5 } },
            { Channel::MuTau, { 4,5 } },
            { Channel::TauTau, { 5,5 } }
        };

        std::cout << boost::format("Processing input file '%1%' into output file '%2%'.\n")
                   % args.input_file() % args.output_file();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto originalTuple = ntuple::CreateEventTuple(args.tree_name(),originalFile.get(),true,ntuple::TreeState::Full);

        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        SummaryInfo* summaryInfo(new SummaryInfo(summaryTuple->data()));
        const Channel channel = Parse<Channel>(args.tree_name());
        const Long64_t n_entries = originalTuple->GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple->GetEntry(current_entry);

            JetOrdering jet_ordering = run_period == Period::Run2017 ? JetOrdering::DeepCSV : JetOrdering::CSV;
            auto event_info =  analysis::MakeEventInfo(channel, originalTuple->data(), run_period, jet_ordering, summaryInfo);
            EventInfoBase& event = *event_info;

            if(event.GetEnergyScale() != EventEnergyScale::Central) continue;
            if(!event.GetTriggerResults().AnyAcceptAndMatch()) continue;
            if(event->extraelec_veto || event->extramuon_veto) continue;

            if(!genMatchInfo.count(channel))
                throw exception("Gen Match numbers not found for this channel.");
            std::pair<Float_t,Float_t> genMatchNumbers = genMatchInfo.at(channel);
            if(event->gen_match_1 == genMatchNumbers.first && event->gen_match_2 == genMatchNumbers.second)
                anaData.mass_bothGood().Fill((event->p4_1 + event->p4_2).mass());
            if(event->gen_match_1 != genMatchNumbers.first && event->gen_match_2 != genMatchNumbers.second)
                anaData.mass_bothBad().Fill((event->p4_1 + event->p4_2).mass());
            if(event->gen_match_1 == genMatchNumbers.first && event->gen_match_2 != genMatchNumbers.second)
                anaData.mass_1stGood_2ndBad().Fill((event->p4_1 + event->p4_2).mass());
            if(event->gen_match_1 != genMatchNumbers.first && event->gen_match_2 == genMatchNumbers.second)
                anaData.mass_1stBad_2ndGood().Fill((event->p4_1 + event->p4_2).mass());
            if(event->gen_match_1 != genMatchNumbers.first || event->gen_match_2 != genMatchNumbers.second)
                anaData.mass_allBad().Fill((event->p4_1 + event->p4_2).mass());

        }

    }

private:
    Arguments args;
    std::shared_ptr<TFile> output;
    GenMatchCheckData anaData;
    analysis::Period run_period;
};

} // namespace analysis

PROGRAM_MAIN(analysis::GenMatchCheck, Arguments)
