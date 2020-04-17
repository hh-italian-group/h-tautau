/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Core/include/ProgressReporter.h"
#include "AnalysisTools/Run/include/program_main.h"
#include "h-tautau/Analysis/include/EventCacheProvider.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Core/include/CacheTuple.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, channels);
    REQ_ARG(analysis::Period, period);
    REQ_ARG(std::string, selections);
    REQ_ARG(std::string, unc_sources);
    REQ_ARG(std::string, btaggers);
    REQ_ARG(bool, isData);
    REQ_ARG(bool, hasBjetPair);
    OPT_ARG(bool, runSVFit, true);
    OPT_ARG(bool, runKinFit, true);
    OPT_ARG(Long64_t, max_events_per_tree, std::numeric_limits<Long64_t>::max());
    OPT_ARG(Long64_t, begin_entry_index, 0);
    OPT_ARG(Long64_t, end_entry_index, std::numeric_limits<Long64_t>::max());
    OPT_ARG(std::string, working_path, "./");
    OPT_ARG(analysis::DiscriminatorWP, btag_wp, analysis::DiscriminatorWP::Medium);
};

namespace analysis {

class CacheTupleProducer {
public:
    using CacheEvent = cache_tuple::CacheEvent;
    using CacheTuple = cache_tuple::CacheTuple;
    using CacheSummaryTuple = cache_ntuple::CacheSummaryTuple;
    using clock = std::chrono::system_clock;

    CacheTupleProducer(const Arguments& _args) :
            args(_args), outputFile(root_ext::CreateRootFile(args.output_file())),
            cacheSummary("summary", outputFile.get(), false), start(clock::now()),
            progressReporter(10, std::cout)
    {
        const auto signal_modes = SplitValueListT<SignalMode>(args.selections(), false, ",");
        for(auto signal_mode : signal_modes)
            signalObjectSelectors.emplace_back(signal_mode);

        EventCandidate::InitializeUncertainties(args.period(), false, args.working_path(),
                                                signalObjectSelectors.at(0).GetTauVSjetDiscriminator().first);

        unc_sources = SplitValueListT<UncertaintySource>(args.unc_sources(), false, ",");
        channels = SplitValueListT<Channel>(args.channels(), false, ",");
        auto btagger_kinds = SplitValueListT<BTaggerKind>(args.btaggers(), false, ",");
        for(BTaggerKind kind : btagger_kinds) {
            if(kind == BTaggerKind::HHbtag)
                throw exception("CacheTupleProducer: HHbtag is not supported.");
            btaggers.emplace_back(args.period(), kind);
        }

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
        std::map<Channel, std::shared_ptr<ntuple::EventTuple>> map_event;
        for(Channel channel : channels){
            try {
                auto originalTuple = ntuple::CreateEventTuple(ToString(channel), originalFile.get(), true,
                                                              ntuple::TreeState::Full);
                const Long64_t n_entries = std::min(originalTuple->GetEntries(), args.end_entry_index())
                                           - args.begin_entry_index();
                const Long64_t n_events = std::min(args.max_events_per_tree(), n_entries);
                map_event[channel] = originalTuple;
                n_tot_events += static_cast<size_t>(n_events);
            } catch(std::runtime_error) {
            }
        }

        size_t n_processed_events = 0;
        progressReporter.SetTotalNumberOfEvents(n_tot_events);
        for(Channel channel : channels){
            if(!map_event.count(channel)) {
                std::cout << "Channel: " << channel << " not found." << std::endl;
                continue;
            }
            std::cout << "Channel: " << channel << std::endl;
            CacheTuple cache(ToString(channel), outputFile.get(), false);
            cache.SetAutoFlush(1000);
            cache.SetMaxVirtualSize(10000000);
            auto& originalTuple = *map_event.at(channel);
            const Long64_t n_entries = originalTuple.GetEntries();
            Long64_t n_processed_events_channel = 0;
            for(Long64_t current_entry = args.begin_entry_index();
                    current_entry < n_entries && n_processed_events_channel < args.max_events_per_tree()
                    && current_entry < args.end_entry_index(); ++current_entry) {
                originalTuple.GetEntry(current_entry);
                originalTuple().isData = args.isData();
                originalTuple().period = static_cast<int>(args.period());
                FillCacheTuple(cache, current_entry, originalTuple.data());
                ++n_processed_events_channel;
                ++n_processed_events;
                if(n_processed_events % 100 == 0) progressReporter.Report(n_processed_events, false);
            }
            progressReporter.Report(n_processed_events, true);
            cache.Write();
            const auto stop = clock::now();
            const auto exeTime = std::chrono::duration_cast<std::chrono::seconds>(stop - start).count();
            cacheSummary().exeTime = static_cast<UInt_t>(exeTime);
            cacheSummary.Fill();
            cacheSummary.Write();
        }
        progressReporter.Report(n_tot_events, true);
    }

private:
    void FillCacheTuple(CacheTuple& cacheTuple, Long64_t original_entry, const ntuple::Event& event)
    {
        if(!SignalObjectSelector::PassLeptonVetoSelection(event)) return;
        if(!SignalObjectSelector::PassMETfilters(event, args.period(), args.isData())) return;

        EventCacheProvider cache_provider;
        for(auto [unc_source, unc_scale] : EnumerateUncVariations(unc_sources)) {
            auto event_candidate = std::make_shared<EventCandidate>(event, unc_source, unc_scale);
            if(unc_scale != UncertaintyScale::Central && event_candidate->IsSameAsCentral()) continue;
            std::set<size_t> htt_indices;
            std::set<std::pair<size_t, size_t>> hh_indices;
            for(const SignalObjectSelector& signalObjectSelector : signalObjectSelectors) {
                for(const BTagger& bTagger : btaggers) {
                    const auto event_info = EventInfo::Create(event_candidate, signalObjectSelector, bTagger,
                                                              args.btag_wp());
                    if(!event_info) continue;
                    if(args.hasBjetPair() && !event_info->HasBjetPair()) continue;

                    const size_t htt_index = event_info->GetHttIndex();
                    if(args.runSVFit() && !htt_indices.count(htt_index)) {
                        cacheSummary().numberOfTimesSVFit++;
                        const sv_fit_ana::FitResults& fit_results = event_info->GetSVFitResults(true);
                        cache_provider.AddSVfitResults(htt_index, unc_source, unc_scale, fit_results);
                        htt_indices.insert(htt_index);
                    }

                    if(args.runKinFit() && event_info->HasBjetPair()) {
                        const size_t hbb_index = event_info->GetSelectedSignalJets().bjet_pair.ToIndex();
                        const auto hh_pair = std::make_pair(htt_index, hbb_index);
                        if(!hh_indices.count(hh_pair)) {
                            cacheSummary().numberOfTimesKinFit++;
                            const kin_fit::FitResults& fit_results = event_info->GetKinFitResults(true);
                            cache_provider.AddKinFitResults(htt_index, hbb_index, unc_source, unc_scale, fit_results);
                            hh_indices.insert(hh_pair);
                        }
                    }
                }
            }
        }

        if(!cache_provider.IsEmpty()) {
            cacheTuple().run = event.run;
            cacheTuple().lumi = event.lumi;
            cacheTuple().evt = event.evt;
            cacheTuple().entry_index = original_entry;
            cache_provider.FillEvent(cacheTuple());
            cacheTuple.Fill();
        }
    }

private:
    Arguments args;
    std::shared_ptr<TFile> outputFile;
    CacheSummaryTuple cacheSummary;
    const clock::time_point start;
    std::vector<Channel> channels;
    std::vector<SignalObjectSelector> signalObjectSelectors;
    std::vector<UncertaintySource> unc_sources;
    std::vector<BTagger> btaggers;
    tools::ProgressReporter progressReporter;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CacheTupleProducer, Arguments)
