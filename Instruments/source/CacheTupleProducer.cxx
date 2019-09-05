/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"

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
};

namespace analysis {

class CacheTupleProducer {
public:
    using CacheEvent = cache_tuple::CacheEvent;
    using CacheTuple = cache_tuple::CacheTuple;

    CacheTupleProducer(const Arguments& _args) : args(_args), run_period(Parse<analysis::Period>(args.period()))
    {
        EventCandidate::InitializeJecUncertainty(run_period);
        std::vector<std::string> selections = SplitValueList(args.selections(),false,",");
        for(unsigned n = 0; n < selections.size(); ++n){
            signalObjectSelectors.emplace_back(Parse<analysis::SignalMode>(selections.at(n)));
        }

        std::vector<std::string> sources = SplitValueList(args.unc_sources(),false,",");
        for(unsigned n = 0; n < sources.size(); ++n){
            unc_sources.emplace_back(Parse<analysis::UncertaintySource>(sources.at(n)));
        }

        std::vector<std::string> jet_orderings = SplitValueList(args.jet_orderings(),false,",");
        for(unsigned n = 0; n < jet_orderings.size(); ++n){
            vector_jet_ordering.emplace_back(Parse<JetOrdering>(jet_orderings.at(n)));
        }

        std::vector<std::string> vector_channel = SplitValueList(args.channels(),false,",");
        for(unsigned n = 0; n < vector_channel.size(); ++n){
            channels.push_back(vector_channel.at(n));
        }

    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% selection.\n")
                   % args.input_file() % args.output_file() % args.selections();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        for(unsigned c = 0; c < channels.size(); ++c){
            auto originalTuple = ntuple::CreateEventTuple(channels.at(c),originalFile.get(),true,ntuple::TreeState::Full);
            CacheTuple cache(channels.at(c), outputFile.get(), false);
            for(const auto& event : *originalTuple) {
                FillCacheTuple(cache, event);
                cache.Fill();
            }
            cache.Write();
        }




    }

private:

    void FillCacheTuple(CacheTuple& cacheTuple, const ntuple::Event& event) const
    {

        std::set<size_t> Htt_indexes;
        std::set<std::pair<size_t,size_t>> HH_indexes;

        for(unsigned n = 0; n < signalObjectSelectors.size(); ++n){
            for(unsigned m = 0; m < vector_jet_ordering.size(); ++m){
                for(unsigned h = 0; h < unc_sources.size(); ++h){
                    for(int variation = -1; variation < 2; ++variation){
                        cacheTuple().run = event.run;
                        cacheTuple().lumi = event.lumi;
                        cacheTuple().evt = event.evt;
                        SignalObjectSelector signalObjectSelector = signalObjectSelectors.at(n);
                        JetOrdering jet_ordering = vector_jet_ordering.at(m);
                        UncertaintySource unc_source = unc_sources.at(h);
                        UncertaintyScale scale = static_cast<UncertaintyScale>(variation);
                        if(scale != UncertaintyScale::Central && unc_source == UncertaintySource::None) continue;
                        if(scale == UncertaintyScale::Central && unc_source != UncertaintySource::None) continue;
                        boost::optional<EventInfoBase> event_info_base = CreateEventInfo(event,signalObjectSelector,nullptr,run_period,jet_ordering,unc_source,scale);
                        if(!event_info_base.is_initialized()) continue;


                        if(args.hasBjetPair() && !event_info_base->HasBjetPair()) continue;
                        if(!signalObjectSelector.PassLeptonVetoSelection(event)) continue;
                        if(!signalObjectSelector.PassMETfilters(event,run_period,args.isData())) continue;

                        size_t selected_htt_index = event_info_base->GetHttIndex();
                        size_t selected_hbb_index = ntuple::CombinationPairToIndex(event_info_base->GetSelectedSignalJets().selectedBjetPair);
                        if(!Htt_indexes.count(selected_htt_index)){
                            const sv_fit_ana::FitResults& result = event_info_base->GetSVFitResults();
                            cacheTuple().SVfit_Higgs_index.push_back(selected_htt_index);
                            cacheTuple().SVfit_is_valid.push_back(result.has_valid_momentum);
                            cacheTuple().SVfit_p4.push_back(result.momentum);
                            cacheTuple().SVfit_p4_error.push_back(LorentzVectorM(result.momentum_error));
                            cacheTuple().SVfit_mt.push_back(static_cast<Float_t>(result.transverseMass));
                            cacheTuple().SVfit_mt_error.push_back(static_cast<Float_t>(result.transverseMass_error));
                            cacheTuple().SVfit_unc_source.push_back(static_cast<Int_t>(unc_sources.at(h)));
                            cacheTuple().SVfit_unc_scale.push_back(variation);
                            Htt_indexes.insert(selected_htt_index);
                        }


                        std::pair<size_t,size_t> hh_pair = std::make_pair(selected_htt_index,selected_hbb_index);
                        if(!HH_indexes.count(hh_pair)){
                            const kin_fit::FitResults& result = event_info_base->GetKinFitResults();
                            cacheTuple().kinFit_Higgs_index.push_back(selected_htt_index);
                            cacheTuple().kinFit_jetPairId.push_back(selected_hbb_index);
                            cacheTuple().kinFit_m.push_back(static_cast<Float_t>(result.mass));
                            cacheTuple().kinFit_chi2.push_back(static_cast<Float_t>(result.chi2));
                            cacheTuple().kinFit_convergence.push_back(result.convergence);
                            cacheTuple().kinFit_unc_source.push_back(static_cast<Int_t>(unc_sources.at(h)));
                            cacheTuple().kinFit_unc_scale.push_back(variation);
                            HH_indexes.insert(hh_pair);
                        }

                    }
                }
            }

        }




    }

private:
    Arguments args;
    analysis::Period run_period;
    std::vector<std::string> channels;
    std::vector<SignalObjectSelector> signalObjectSelectors;
    std::vector<UncertaintySource> unc_sources;
    std::vector<JetOrdering> vector_jet_ordering;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CacheTupleProducer, Arguments)
