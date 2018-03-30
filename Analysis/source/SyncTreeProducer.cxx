/*! Produce synchronization tree.
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

struct Arguments {
    REQ_ARG(std::string, mode);
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, period);
    REQ_ARG(std::string, output_file);
    OPT_ARG(std::string, sample_type, "signal");
};

namespace analysis {

enum class SyncMode { HTT, HH };

ENUM_NAMES(SyncMode) = {
    { SyncMode::HTT, "htt" },
    { SyncMode::HH, "hh" }
};

class SyncTreeProducer {
public:
    using Event = ntuple::Event;
    using EventTuple = ntuple::EventTuple;
    using SyncEvent = htt_sync::SyncEvent;
    using SyncTuple = htt_sync::SyncTuple;

    static constexpr float default_value = std::numeric_limits<float>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    SyncTreeProducer(const Arguments& _args) : args(_args), eventWeights(Period::Run2016, DiscriminatorWP::Medium)
    {
        std::istringstream ss_mode(args.mode());
        ss_mode >> syncMode;
        run_period = analysis::EnumNameMap<analysis::Period>::GetDefault().Parse(args.period());
    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% mode.\n")
                   % args.input_file() % args.output_file() % args.mode();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        auto originalTuple = ntuple::CreateEventTuple(args.tree_name(),originalFile.get(),true,ntuple::TreeState::Full);
        SyncTuple sync(args.tree_name(), outputFile.get(), false);
        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        std::shared_ptr<SummaryInfo> summaryInfo(new SummaryInfo(summaryTuple->data()));
        const Channel channel = Parse<Channel>(args.tree_name());
        const Long64_t n_entries = originalTuple->GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple->GetEntry(current_entry);

            ntuple::JetPair bjet_pair;
            if(run_period == Period::Run2016)
                bjet_pair = EventInfoBase::SelectBjetPair(originalTuple->data(), cuts::btag_2016::pt,
                                                                 cuts::btag_2016::eta, JetOrdering::CSV);
            if(run_period == Period::Run2017)
                bjet_pair = EventInfoBase::SelectBjetPair(originalTuple->data(), cuts::btag_2017::pt,
                                                          cuts::btag_2017::eta, JetOrdering::DeepCSV);
            auto eventInfoPtr = MakeEventInfo(channel, originalTuple->data(), bjet_pair, summaryInfo.get());
            EventInfoBase& event = *eventInfoPtr;

            if(event.GetEnergyScale() != EventEnergyScale::Central) continue;

            if(syncMode == SyncMode::HH) {
                if(event->extraelec_veto || event->extramuon_veto) continue;
            }

            htt_sync::FillSyncTuple(event,sync,run_period);
        }

        sync.Write();
    }

private:
    Arguments args;
    SyncMode syncMode;
    analysis::Period run_period;
    mc_corrections::EventWeights eventWeights;
};

} // namespace analysis

PROGRAM_MAIN(analysis::SyncTreeProducer, Arguments)
