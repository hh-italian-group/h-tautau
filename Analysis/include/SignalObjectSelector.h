/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2017.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "h-tautau/Cuts/include/H_tautau_2017_baseline.h"
#include "h-tautau/Analysis/include/MetFilters.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Analysis/include/EventCandidate.h"

namespace analysis {

enum class SignalMode { HTT = 1, HTT_sync = 2, TauPOG_default = 3, TauPOG_deepTauVsJet = 4,
                        TauPOG_deepTauVsJet_full = 5, TauPOG_dpfTau = 6, HH_legacy = 7, HH = 8,
                        Skimmer = 9, TauPOG_Skimmer = 10 };

ENUM_NAMES(SignalMode) = {
    { SignalMode::HTT, "HTT" },
    { SignalMode::HTT_sync, "HTT_sync" },
    { SignalMode::TauPOG_default, "TauPOG_default" },
    { SignalMode::TauPOG_deepTauVsJet, "TauPOG_deepTauVsJet" },
    { SignalMode::TauPOG_deepTauVsJet_full, "TauPOG_deepTauVsJet_full" },
    { SignalMode::HH_legacy, "HH_legacy" },
    { SignalMode::HH, "HH" },
    { SignalMode::Skimmer, "Skimmer" },
    { SignalMode::TauPOG_Skimmer, "TauPOG_Skimmer" }
};

namespace jet_ordering {

    template<typename LorentzVector>
    struct JetInfo {
        LorentzVector p4;
        size_t index;
        double tag;

        JetInfo() : index(0), tag(0.0) { }

        JetInfo(const LorentzVector _p4, size_t _index, double _tag)
            : p4(_p4), index(_index), tag(_tag) { }

    };

    template<typename LorentzVector>
    bool CompareJets(const JetInfo<LorentzVector>& jet_1, const JetInfo<LorentzVector>& jet_2,
                     double pt_thr, double eta_thr)
    {
        const auto eta1 = std::abs(jet_1.p4.eta());
        const auto eta2 = std::abs(jet_2.p4.eta());
        if(eta1 < eta_thr && eta2 >= eta_thr) return true;
        if(eta1 >= eta_thr && eta2 < eta_thr) return false;
        const auto pt1 = jet_1.p4.pt();
        const auto pt2 = jet_2.p4.pt();
        if(pt1 > pt_thr && pt2 <= pt_thr) return true;
        if(pt1 <= pt_thr && pt2 > pt_thr) return false;

        if(jet_1.tag != jet_2.tag)
            return jet_1.tag > jet_2.tag;
        return pt1 > pt2;
    };

    template<typename LorentzVector>
    std::vector<JetInfo<LorentzVector>> OrderJets(
                                  const std::vector<JetInfo<LorentzVector>>& jet_info_vector,
                                  bool apply_hard_cut,
                                  double pt_cut = std::numeric_limits<double>::lowest(),
                                  double eta_cut = std::numeric_limits<double>::max())
    {
        const auto comparitor = [&](const JetInfo<LorentzVector>& jet_1,
                                    const JetInfo<LorentzVector>& jet_2) -> bool {
            return analysis::jet_ordering::CompareJets(jet_1, jet_2, pt_cut, eta_cut);
        };

        std::vector<JetInfo<LorentzVector>> jets_ordered;
        if(apply_hard_cut){
            for(size_t n = 0; n < jet_info_vector.size(); ++n) {
                if(jet_info_vector.at(n).p4.Pt() > pt_cut && std::abs(jet_info_vector.at(n).p4.eta()) < eta_cut)
                    jets_ordered.push_back(jet_info_vector.at(n));
            }
        }
        else
            jets_ordered = jet_info_vector;
        std::sort(jets_ordered.begin(), jets_ordered.end(), comparitor);
        return jets_ordered;
    }

}


class SignalObjectSelector {
public:
    using LegPair = ntuple::LegPair;
    using LepCandidate = LeptonCandidate<ntuple::TupleLepton>;

    SignalObjectSelector(SignalMode _mode);

    bool PassLeptonSelection(const LepCandidate& lepton, Channel channel, const size_t legId) const;
    boost::optional<size_t> GetHiggsCandidateIndex(EventCandidate& event_candidate) const;
    bool PassLeptonVetoSelection(const ntuple::Event& event) const;
    bool PassMETfilters(const ntuple::Event& event, Period period, bool is_Data) const;

    struct SelectedSignalJets{
        LegPair selectedBjetPair;
        LegPair selectedVBFjetPair;
        size_t n_bjets;

        SelectedSignalJets();
        bool HasBjetPair(size_t njets) const;
        bool HasVBFPair(size_t njets) const;
        bool isSelectedBjet(size_t n) const;
        bool isSelectedVBFjet(size_t n) const;
    };

    static SelectedSignalJets SelectSignalJets(EventCandidate& event_candidate,
                                               const analysis::Period& period,
                                               analysis::JetOrdering jet_ordering,
                                               size_t selected_higgs_index,
                                               analysis::UncertaintySource unc_source,
                                               analysis::UncertaintyScale unc_scale);

     template<typename LorentzVector>
     static bool PassEcalNoiceVetoJets(const LorentzVector& jet_p4, Period period, DiscriminatorIdResults jets_pu_id)
     {
         if(period !=  analysis::Period::Run2017)
             return true;

         const double abs_eta = std::abs(jet_p4.eta());
         return !(jet_p4.pt() < cuts::hh_bbtautau_2017::jetID::max_pt_veto &&
                     abs_eta > cuts::hh_bbtautau_2017::jetID::eta_low_veto &&
                     abs_eta < cuts::hh_bbtautau_2017::jetID::eta_high_veto && !jets_pu_id.Passed(analysis::DiscriminatorWP::Loose));
     }

private:
    bool PassHTT_LeptonSelection(const LepCandidate& lepton, Channel channel, bool is_sync) const;
    bool PassTauPOG_LeptonSelection(const LepCandidate& lepton, Channel channel) const;
    bool PassHH_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId) const;
    bool PassHH_legacy_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId) const;
    bool PassSkimmer_LeptonSelection(const LepCandidate& lepton) const;
    bool PassTauPOG_Skimmer_LeptonSelection(const LepCandidate& lepton) const;



private:
    SignalMode mode;
    double DR2_leptons;
    TauIdDiscriminator discriminator;
};

} // namespace analysis
