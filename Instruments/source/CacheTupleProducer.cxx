/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/McCorrections/include/EventWeights.h"
#include "hh-bbtautau/Analysis/include/SampleDescriptorConfigEntryReader.h"
#include "h-tautau/Core/include/CacheTuple.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Analysis/include/SignalObjectSelector.h"

struct Arguments {
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, output_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, trigger_cfg);
    REQ_ARG(std::string, period);
    REQ_ARG(std::string, selection);
    REQ_ARG(std::string, unc_source);
    REQ_ARG(bool, isData);
};

namespace analysis {

class CacheTupleProducer {
public:
    using CacheEvent = cache_tuple::CacheEvent;
    using CacheTuple = cache_tuple::CacheTuple;

    CacheTupleProducer(const Arguments& _args) : args(_args), run_period(Parse<analysis::Period>(args.period()))
    {
        std::vector<std::string> selections = SplitValueList(args.selection(),false,",");
        for(unsigned n = 0; n < selections.size(); ++n){
            SignalObjectSelector signalObjectSelector(Parse<analysis::SignalMode>(selections.at(n)));
            signalObjectSelectors.push_back(signalObjectSelector);
        }

        std::vector<std::string> sources = SplitValueList(args.unc_source(),false,",");
        for(unsigned n = 0; n < sources.size(); ++n){
            UncertaintySource source = Parse<analysis::UncertaintySource>(sources.at(n));
            unc_sources.push_back(source);
        }
    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% selection.\n")
                   % args.input_file() % args.output_file() % args.selection();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        auto originalTuple = ntuple::CreateEventTuple(args.tree_name(),originalFile.get(),true,ntuple::TreeState::Full);
        CacheTuple cache(args.tree_name(), outputFile.get(), false);
        auto summaryTuple = ntuple::CreateSummaryTuple("summary", originalFile.get(), true, ntuple::TreeState::Full);
        summaryTuple->GetEntry(0);
        SummaryInfo summaryInfo(summaryTuple->data(), Parse<Channel>(args.tree_name()), args.trigger_cfg());
        for(const auto& event : *originalTuple) {
            FillCacheTuple(cache, event,summaryInfo);
            cache.Fill();
        }

        cache.Write();

    }

private:

    void FillCacheTuple(CacheTuple& cacheTuple, const ntuple::Event& event,const SummaryInfo& summaryInfo) const
    {
        static const std::map<Channel, std::vector<std::string>> triggerPaths = {
            { Channel::ETau, { "HLT_Ele32_WPTight_Gsf_v", "HLT_Ele35_WPTight_Gsf_v", "HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v" } },
            { Channel::MuTau, { "HLT_IsoMu24_v", "HLT_IsoMu27_v", "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v" } },
            { Channel::TauTau, { "HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v",
                "HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v",
                "HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"} },
            { Channel::MuMu, { "HLT_IsoMu24_v", "HLT_IsoMu27_v" } },
        };
        const Channel channel = Parse<Channel>(args.tree_name());

        if(event.extraelec_veto || event.extramuon_veto) return;

        //JetOrdering jet_ordering = run_period == Period::Run2017 ? JetOrdering::DeepCSV : JetOrdering::CSV;
        JetOrdering jet_ordering = JetOrdering::DeepFlavour;
        std::set<size_t> Htt_indexes;
        std::set<std::pair<size_t,size_t>> HH_indexes;
        for(unsigned n = 0; n < signalObjectSelectors.size(); ++n){
            for(unsigned h = 0; h < unc_sources.size(); ++h){
                for(int l = -1; l < 2; ++l){
                    SignalObjectSelector signalObjectSelector = signalObjectSelectors.at(n);
                    boost::optional<EventInfoBase> event_info_base = CreateEventInfo(event,signalObjectSelector,&summaryInfo,run_period,jet_ordering);
                    if(!event_info_base.is_initialized()) continue;
                    if(!event_info_base->GetTriggerResults().AnyAcceptAndMatchEx(triggerPaths.at(channel), event_info_base->GetFirstLeg().GetMomentum().pt(),
                                                                                                        event_info_base->GetSecondLeg().GetMomentum().pt())) continue;
                    if(!event_info_base->HasBjetPair()) continue;
                    if(!signalObjectSelector.PassLeptonVetoSelection(event)) continue;
                    if(!signalObjectSelector.PassMETfilters(event,run_period,args.isData())) continue;
                    for(size_t leg_id = 1; leg_id <= 2; ++leg_id) {
                        const LepCandidate& lepton = event_info_base->GetLeg(leg_id);
                        if(lepton->leg_type() == LegType::tau){
                            if(!lepton->Passed(TauIdDiscriminator::byDeepTau2017v2p1VSjet, DiscriminatorWP::Medium)) continue;
                        }
                    }

                    size_t selected_htt_index = event_info_base->GetHttIndex();
                    size_t selected_hbb_index = ntuple::CombinationPairToIndex(event_info_base->GetSelectedSignalJets().selectedBjetPair);
                    if(!Htt_indexes.count(selected_htt_index)){
                        const sv_fit_ana::FitResults& result = event_info_base->GetSVFitResults();
                        cacheTuple().SVfit_Higgs_index.push_back(selected_htt_index);
                        cacheTuple().SVfit_is_valid.push_back(result.has_valid_momentum);
                        cacheTuple().SVfit_p4.push_back(LorentzVectorM(result.momentum));
                        cacheTuple().SVfit_p4_error.push_back(LorentzVectorM(result.momentum_error));
                        cacheTuple().SVfit_mt.push_back(result.transverseMass);
                        cacheTuple().SVfit_mt_error.push_back(result.transverseMass_error);
                        cacheTuple().SVfit_unc_source.push_back(static_cast<Int_t>(unc_sources.at(h)));
                        cacheTuple().SVfit_unc_scale.push_back(l);
                        Htt_indexes.insert(selected_htt_index);
                    }


                    std::pair<size_t,size_t> hh_pair = std::make_pair(selected_htt_index,selected_hbb_index);
                    if(!HH_indexes.count(hh_pair)){
                        const kin_fit::FitResults& result = event_info_base->GetKinFitResults();
                        cacheTuple().kinFit_Higgs_index.push_back(selected_htt_index);
                        cacheTuple().kinFit_jetPairId.push_back(selected_hbb_index);
                        cacheTuple().kinFit_m.push_back(result.mass);
                        cacheTuple().kinFit_chi2.push_back(result.chi2);
                        cacheTuple().kinFit_convergence.push_back(result.convergence);
                        cacheTuple().kinFit_unc_source.push_back(static_cast<Int_t>(unc_sources.at(h)));
                        cacheTuple().kinFit_unc_scale.push_back(l);
                        HH_indexes.insert(hh_pair);
                    }

                }
            }
        }




    }

private:
    Arguments args;
    analysis::Period run_period;
    std::vector<SignalObjectSelector> signalObjectSelectors;
    std::vector<UncertaintySource> unc_sources;
};

} // namespace analysis

PROGRAM_MAIN(analysis::CacheTupleProducer, Arguments)
