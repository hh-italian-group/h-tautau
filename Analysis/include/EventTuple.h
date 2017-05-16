/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace ntuple {
using LorentzVectorE = analysis::LorentzVectorE_Float;
using LorentzVectorM = analysis::LorentzVectorM_Float;
using MetCovMatrix = analysis::SquareMatrix<2>;
}


#define LVAR(type, name, n) VAR(type, name##_##n)

#define LEG_DATA(n) \
    LVAR(LorentzVectorM, p4, n) /* 4-momentum */ \
    LVAR(Int_t, q, n) /* Charge */ \
    LVAR(Float_t, dxy, n) /* dxy with respect to primary vertex */ \
    LVAR(Float_t, dz, n) /* dz with respect to primary vertex */ \
    LVAR(Float_t, iso, n) /* MVA iso for hadronic Tau, Delta Beta for muon and electron */ \
    LVAR(Int_t, gen_match, n) /* Generator matching, see Htautau Twiki*/\
    LVAR(LorentzVectorM, gen_p4, n) /* 4-momentum of the matched gen particle */ \
    LVAR(std::vector<uint32_t>, tauId_keys, n) /* keys for tau ID variables */ \
    LVAR(std::vector<float>, tauId_values, n) /* values of tau ID variables */ \
    /**/

#define JVAR(type, name, col) VAR(std::vector<type>, col##_##name)

#define JET_COMMON_DATA(col) \
    JVAR(LorentzVectorE, p4, col) /* Jet 4-momentum */ \
    JVAR(Float_t, csv, col) /* Jet CSV value */ \
    /**/

#define JET_DATA(col) \
    JET_COMMON_DATA(col) \
    JVAR(Float_t, rawf, col) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    JVAR(Float_t, mva, col) /* Jet MVA id value */ \
    JVAR(Int_t, partonFlavour, col) \
    JVAR(Int_t, hadronFlavour, col) \
    /**/

#define FATJET_DATA(col) \
    JET_COMMON_DATA(col) \
    JVAR(Float_t, m_pruned, col) \
    JVAR(Float_t, m_softDrop, col) \
    JVAR(Float_t, n_subjettiness_tau1, col) \
    JVAR(Float_t, n_subjettiness_tau2, col) \
    JVAR(Float_t, n_subjettiness_tau3, col) \
    /**/

#define SUBJET_DATA(col) \
    JET_COMMON_DATA(col) \
    JVAR(size_t, parentIndex, col) \
    /**/

#define MVAR(type, name, col) VAR(type, col##_##name)

#define MET_DATA(col) \
    MVAR(LorentzVectorM, p4, col) /* MET 4-momentum */ \
    MVAR(MetCovMatrix, cov, col) /* pf met covariance matrix */ \
    /**/

#define EVENT_DATA() \
    VAR(UInt_t, run) /* run */ \
    VAR(UInt_t, lumi) /* lumi section */ \
    VAR(ULong64_t, evt) /* event number */ \
    VAR(Int_t, channelId) /* Channel: eTau, muTau or tauTau */ \
    VAR(Int_t, eventEnergyScale) /* event type category */ \
    VAR(Int_t, genEventType) /* gen event type */ \
    VAR(Float_t, genEventWeight) /* gen event weight */ \
    VAR(UInt_t, storageMode) /* for non-central ES, description of the relation with central ES event */ \
	/* Event Weights Variables */ \
    VAR(Float_t, weight_btag) \
    VAR(Float_t, weight_lepton) \
    VAR(Float_t, weight_ttbar_pt) \
    VAR(Float_t, weight_ttbar_merge) \
    VAR(Float_t, weight_PU) \
    VAR(Float_t, weight_dy) \
    VAR(Float_t, weight_sm) \
    /* VAR(Float_t, shape_denominator_weight) */ \
    /* Event Variables */ \
    VAR(Int_t, npv) /* NPV */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
	VAR(UInt_t, n_jets) /* Number of jets in the event */\
    VAR(Float_t, ht_other_jets) /* Ht of all jets in the event except the first 2 jets */\
    /* Trigger results */ \
    VAR(ULong64_t, trigger_accepts) /* Trigger accept bits for the selected triggers */ \
    VAR(ULong64_t, trigger_matches) /* Leg matching results for the selected triggers */ \
    /* SV Fit variables */ \
    VAR(LorentzVectorM, SVfit_p4) /* SVfit using integration method */ \
    VAR(Float_t, SVfit_mt) /* SVfit using integration method */ \
    /* Signal leptons */ \
    LEG_DATA(1) /* muon for muTau, electron for eTau, electron for eMu, Leading (in pT) Tau for tauTau */ \
    LEG_DATA(2) /* hadronic Tau for muTau and eTau, Muon for eMu, Trailing (in pT) Tau for tauTau */ \
    /* Met related variables */ \
    MET_DATA(pfMET) \
    VAR(UInt_t, metFilters) \
    /* Candidate Jets: jets after applying Jet energy corrections (excluding hadronic Tau) */ \
    JET_DATA(jets) \
    FATJET_DATA(fatJets) \
    SUBJET_DATA(subJets) \
    /* KinFit Variables */ \
    VAR(std::vector<UInt_t>, kinFit_jetPairId) /* indices of jet pairs for which KinFit is calculated */\
    VAR(std::vector<Float_t>, kinFit_m) /* KinFit m_bbtt mass */\
    VAR(std::vector<Float_t>, kinFit_chi2) /*  KinFit chi2 value*/ \
    VAR(std::vector<Int_t>, kinFit_convergence) /* KinFit convergence code */\
    /* Generator level information */\
    VAR(UInt_t, lhe_n_partons) \
    VAR(UInt_t, lhe_n_c_partons) \
    VAR(UInt_t, lhe_n_b_partons) \
    VAR(Float_t, lhe_HT) \
    VAR(Float_t, lhe_H_m) \
    VAR(Float_t, lhe_hh_m) \
    VAR(Float_t, lhe_hh_cosTheta) \
    VAR(std::vector<Int_t>, genParticles_pdg) \
    VAR(std::vector<LorentzVectorM>, genParticles_p4) \
    VAR(UInt_t, genParticles_nPromptElectrons) \
    VAR(UInt_t, genParticles_nPromptMuons) \
    VAR(UInt_t, genParticles_nPromptTaus) \
    VAR(UInt_t, genJets_nTotal) \
    VAR(UInt_t, jets_nTotal_partonFlavour_b) \
    VAR(UInt_t, jets_nTotal_partonFlavour_c) \
    VAR(UInt_t, jets_nTotal_hadronFlavour_b) \
    VAR(UInt_t, jets_nTotal_hadronFlavour_c) \
    VAR(std::vector<LorentzVectorE>, genJets_p4) \
    VAR(std::vector<Int_t>, genJets_partonFlavour) \
    VAR(std::vector<Int_t>, genJets_hadronFlavour) \
    /* Vetos */\
    VAR(Bool_t, trigger_match) /* True if event passes trigger match. */ \
    VAR(Bool_t, dilepton_veto) /* Event is vetoed by the dilepton veto if true */ \
    VAR(Bool_t, extraelec_veto) /* Event is vetoed by the extra electron veto if true */ \
    VAR(Bool_t, extramuon_veto) /* Event is vetoed by the extra muon veto if true */ \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(ntuple, Event, EventTuple, EVENT_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(ntuple, EventTuple, EVENT_DATA)
#undef VAR
#undef EVENT_DATA
#undef LEG_DATA
#undef LVAR
#undef JET_DATA
#undef FATJET_DATA
#undef SUBJET_DATA
#undef JET_COMMON_DATA
#undef JVAR
#undef MET_DATA
#undef MVAR

namespace ntuple {
template<typename T>
constexpr T DefaultFillValue() { return std::numeric_limits<T>::lowest(); }

using JetPair = std::pair<size_t, size_t>;

inline size_t NumberOfCombinationPairs(size_t n_jets) { return n_jets * (n_jets - 1); }

inline size_t CombinationPairToIndex(const JetPair& pair, size_t n_jets)
{
    const size_t min = std::min(pair.first, pair.second);
    const size_t max = std::max(pair.first, pair.second);
    if(n_jets < 2 || min == max || max >= n_jets)
        throw analysis::exception("bad combination pair (%1%, %2%) for n b-jets = %3%.")
            % pair.first % pair.second % n_jets;
    const size_t corr = pair.first < pair.second ? -1 : 0;
    return pair.first * (n_jets - 1) + pair.second + corr;
}

inline JetPair CombinationIndexToPair(size_t index, size_t n_jets)
{
    if(n_jets < 2 || index >= NumberOfCombinationPairs(n_jets))
        throw analysis::exception("bad combination index = %1% for n b-jets = %2%.") % index % n_jets;

    JetPair pair;
    pair.second = index % (n_jets - 1);
    pair.first = (index - pair.second) / (n_jets - 1);
    if(pair.first <= pair.second)
        ++pair.second;
    return pair;
}

inline JetPair UndefinedJetPair()
{
    static JetPair pair(std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max());
    return pair;
}


} // namespace ntuple
