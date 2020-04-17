/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"
#include "h-tautau/Analysis/include/MetFilters.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Core/include/Candidate.h"
#include "h-tautau/Analysis/include/EventCandidate.h"

namespace analysis {

enum class SignalMode { HTT, TauPOG, HH, HH_legacy };

ENUM_NAMES(SignalMode) = {
    { SignalMode::HTT, "HTT" },
    { SignalMode::TauPOG, "TauPOG" },
    { SignalMode::HH, "HH" },
    { SignalMode::HH_legacy, "HH_legacy" },
};

namespace jet_ordering {

    struct JetInfo {
        LorentzVectorE p4;
        size_t index;
        double tag;

        JetInfo() : index(0), tag(0.0) {}
        template<typename LVector>
        JetInfo(const LVector _p4, size_t _index, double _tag) : p4(_p4), index(_index), tag(_tag) {}
    };

    bool CompareJets(const JetInfo& jet_1, const JetInfo& jet_2, const Cut1D& pt_cut, const Cut1D& eta_cut);
    std::vector<JetInfo> FilterJets(const std::vector<JetInfo>& jets, const Cut1D& pt_cut, const Cut1D& eta_cut);
    std::vector<JetInfo> OrderJets(const std::vector<JetInfo>& jets, bool apply_hard_cut,
                                   const Cut1D& pt_cut, const Cut1D& eta_cut);
}


class SignalObjectSelector {
public:
    using LegPair = ntuple::LegPair;
    using LepCandidate = LeptonCandidate<ntuple::TupleLepton>;
    using JetInfoCollection = std::vector<jet_ordering::JetInfo>;

    struct SelectedSignalJets{
        LegPair bjet_pair;
        LegPair vbf_pair;
        size_t n_bjets;

        SelectedSignalJets();
        bool HasBjetPair() const;
        bool HasVBFPair() const;
        bool isSelectedBjet(size_t n) const;
        bool isSelectedVBFjet(size_t n) const;
    };

    SignalObjectSelector(SignalMode _mode);

    SignalMode GetSignalMode() const { return mode; }

    bool PassLeptonSelection(const LepCandidate& lepton, Channel channel, const size_t legId,
                             bool is_sync = false) const;
    boost::optional<size_t> GetHiggsCandidateIndex(const EventCandidate& event_candidate, bool is_sync = false) const;
    static bool PassLeptonVetoSelection(const ntuple::Event& event);
    static bool PassMETfilters(const ntuple::Event& event, Period period, bool is_Data);
    std::pair<TauIdDiscriminator, DiscriminatorWP> GetTauVSjetDiscriminator() const;
    std::pair<TauIdDiscriminator, DiscriminatorWP> GetTauVSeDiscriminator(Channel channel) const;
    std::pair<TauIdDiscriminator, DiscriminatorWP> GetTauVSmuDiscriminator(Channel channel) const;
    std::pair<DiscriminatorWP, DiscriminatorWP>  GetTauVSjetSidebandWPRange() const;

    static JetInfoCollection CreateJetInfos(const EventCandidate& event_candidate, const BTagger& bTagger,
                                            bool apply_jet_up_id);
    static JetInfoCollection CreateJetInfos(const EventCandidate& event_candidate, const BTagger& bTagger,
                                            bool apply_jet_up_id, const boost::optional<size_t>& selected_htt_index,
                                            const SelectedSignalJets& selected_signal_jets);
    static SelectedSignalJets SelectSignalJets(const EventCandidate& event_candidate, size_t selected_htt_index,
                                               const BTagger& bTagger, DiscriminatorWP btag_wp);
    static const FatJetCandidate* SelectFatJet(const EventCandidate& event_candidate,
                                               const SelectedSignalJets& selected_signal_jets);

    template<typename LVector>
    static bool PassEcalNoiceVeto(const LVector& jet_p4, Period period, DiscriminatorIdResults jet_pu_id)
    {
        return PassEcalNoiceVetoImpl(LorentzVector(jet_p4), period, jet_pu_id);
    }

private:
    static bool PassEcalNoiceVetoImpl(const LorentzVector& jet_p4, Period period, DiscriminatorIdResults jet_pu_id);

    bool PassHTT_LeptonSelection(const LepCandidate& lepton, Channel channel, bool is_sync = false) const;
    bool PassTauPOG_LeptonSelection(const LepCandidate& lepton, Channel channel) const;
    bool PassHH_LeptonSelection(const LepCandidate& lepton, Channel channel, size_t legId,  bool is_sync = false) const;

private:
    SignalMode mode;
    double DR2_leptons;
};

} // namespace analysis
