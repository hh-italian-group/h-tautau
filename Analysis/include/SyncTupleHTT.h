/*! Definition of a tuple with all event information that is required at the analysis level.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/SmartTree.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "EventInfo.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "AnalysisTypes.h"

#define LVAR(type, name, pref) VAR(type, name##_##pref)
#define JVAR(type, name, suff, pref) VAR(type, suff##name##_##pref)

#define LEG_DATA(pref) \
    LVAR(Float_t, pt, pref) \
    LVAR(Float_t, pt_tau_ES_up, pref) \
    LVAR(Float_t, pt_tau_ES_down, pref) \
    LVAR(Float_t, phi, pref) \
    LVAR(Float_t, eta, pref) \
    LVAR(Float_t, m, pref) \
    LVAR(Float_t, q, pref) \
    LVAR(Float_t, d0, pref) /* is dxy between leading track and first PV */ \
    LVAR(Float_t, dZ, pref) /* dZ between leading track and first PV */ \
    LVAR(Float_t, mt, pref) /* Use MVAMet sqrt(2 * l_pt * met_pt * (1 - cos( d_phi(l, met) )) */ \
    LVAR(Float_t, pfmt, pref) /* As above but using PF Met */ \
    LVAR(Float_t, puppimt, pref) /* As above but using Puppi Met */ \
    LVAR(Float_t, iso, pref) /* iso */ \
    LVAR(Float_t, id_e_mva_nt_loose, pref) /* Non-triggering electron ID MVA score */ \
    LVAR(Float_t, gen_match, pref) /* Type of gen particle matched to object (see HiggsToTauTauWorking2016#MC_Matching) */ \
    LVAR(Float_t, againstElectronLooseMVA6, pref) \
    LVAR(Float_t, againstElectronMediumMVA6, pref) \
    LVAR(Float_t, againstElectronTightMVA6, pref) \
    LVAR(Float_t, againstElectronVLooseMVA6, pref) \
    LVAR(Float_t, againstElectronVTightMVA6, pref) \
    LVAR(Float_t, againstMuonLoose3, pref) \
    LVAR(Float_t, againstMuonTight3, pref) \
    LVAR(Float_t, byCombinedIsolationDeltaBetaCorrRaw3Hits, pref) \
    LVAR(Float_t, byIsolationMVArun2v1DBoldDMwLTraw, pref) \
    LVAR(Float_t, byIsolationMVArun2017v2DBoldDMwLTraw2017, pref) \
    LVAR(Float_t, chargedIsoPtSum, pref) \
    LVAR(Float_t, decayModeFindingOldDMs, pref) \
    LVAR(Float_t, neutralIsoPtSum, pref) \
    LVAR(Float_t, puCorrPtSum, pref) \
    LVAR(Float_t, trigweight, pref) \
    LVAR(Float_t, idisoweight, pref) \
    /**/

#define JET_DATA(suff, pref) \
    JVAR(Float_t, pt, suff, pref) \
    JVAR(Float_t, eta, suff, pref) \
    JVAR(Float_t, phi, suff, pref) \
    JVAR(Float_t, rawf, suff, pref) /* factor to be applied to the jet p4 to obtain its uncorrected p4 */ \
    JVAR(Float_t, csv, suff, pref) \
    JVAR(Float_t, deepcsv, suff, pref) \
    JVAR(Float_t, resolution, suff, pref) /* Jet energy resolution in percentage */ \
    JVAR(Float_t, pt_tau_ES_up, suff, pref) \
    JVAR(Float_t, pt_tau_ES_down, suff, pref) \
    JVAR(Float_t, pt_jet_ES_up, suff, pref) \
    JVAR(Float_t, pt_jet_ES_down, suff, pref) \
    /**/

#define SYNC_DATA() \
    VAR(UInt_t, run) /* Run */ \
    VAR(UInt_t, lumi) /* Lumi */ \
    VAR(ULong64_t, evt) /* Evt */ \
    VAR(UInt_t, sampleId) /* sample id */ \
    /* Event Variables */ \
    VAR(Int_t, npv) /* Number of offline primary vertices */ \
    VAR(Float_t, npu) /* Number of in-time pu interactions added to the event */ \
    VAR(Float_t, rho) /* Use fixedGridRhoFastjetAll */ \
    /* Signal leptons */ \
    LEG_DATA(1) /* Leg 1 (leading tau for tt, electon for et,em muon for mt) */ \
    LEG_DATA(2) /* Leg 2 (trailing tau for tt, tau for et, mt, muon for em) */ \
    /* di-tau system */ \
    VAR(Float_t, pt_tt) /* like HIG-13-004 (p4_l1+p4_l2+p4_MET)_T, use PF met */  \
    VAR(Float_t, mt_tot) /* Use MVA MET (see HiggsToTauTauWorking2016#Synchronisation_Ntuple) */  \
    VAR(Float_t, m_vis) \
    VAR(Float_t, m_sv) /* using MarkovChain MC integration svfit.mass() */ \
    VAR(Float_t, m_sv_tau_ES_up) /*  */ \
    VAR(Float_t, m_sv_tau_ES_down) /*  */ \
    VAR(Float_t, m_sv_jet_ES_up) /*  */ \
    VAR(Float_t, m_sv_jet_ES_down) /*  */ \
    VAR(Float_t, mt_sv) /* using MarkovChain MC integration svfit.transverseMass() */ \
    /* MT2 */ \
    VAR(Double_t, mt2) /* MT2 */ \
    VAR(Double_t, mt2_tau_ES_up) /*  */ \
    VAR(Double_t, mt2_tau_ES_down) /*  */ \
    VAR(Double_t, mt2_jet_ES_up) /*  */ \
    VAR(Double_t, mt2_jet_ES_down) /*  */ \
    /* MET */ \
    VAR(Float_t, met) /* type 1 corrected PF MET */  \
    VAR(Float_t, met_tau_ES_up) /* */  \
    VAR(Float_t, met_tau_ES_down) /* */  \
    VAR(Float_t, met_jet_ES_up) /* */  \
    VAR(Float_t, met_jet_ES_down) /* */  \
    VAR(Float_t, metphi) /* type 1 corrected PF MET phi */  \
    VAR(Float_t, metphi_tau_ES_up) /* */  \
    VAR(Float_t, metphi_tau_ES_down) /* */  \
    VAR(Float_t, metphi_jet_ES_up) /* */  \
    VAR(Float_t, metphi_jet_ES_down) /* */  \
    VAR(Float_t, puppimet) /* Puppi corrrected MET */  \
    VAR(Float_t, puppimetphi) /* Puppi corrected MET phi */  \
    VAR(Float_t, mvamet) /* mva met */  \
    VAR(Float_t, mvametphi) /* mva met phi */  \
    VAR(Float_t, pzetavis) /* see HiggsToTauTauWorking2016#Synchronisation_Ntuple */  \
    VAR(Float_t, pzetamiss) /* use MVA met, see HiggsToTauTauWorking2016#Synchronisation_Ntuple */  \
    VAR(Float_t, pfpzetamiss) /* As above but using pf met */  \
    VAR(Float_t, puppipzetamiss) /* As above but using puppi met */  \
    VAR(Float_t, mvacov00) /* mva met */  \
    VAR(Float_t, mvacov01) /* mva met */  \
    VAR(Float_t, mvacov10) /* mva met */  \
    VAR(Float_t, mvacov11) /* mva met */  \
    VAR(Float_t, metcov00) /* pf met */  \
    VAR(Float_t, metcov01) /* pf met */  \
    VAR(Float_t, metcov10) /* pf met */  \
    VAR(Float_t, metcov11) /* pf met */  \
    /* VBF system (Only fill if njetspt20>=2) */ \
    VAR(Float_t, mjj) /* (jet_1->vector()+jet_2->vector() ).M() */ \
    VAR(Float_t, jdeta) /* delta eta between leading jet and subleading jet */ \
    VAR(Float_t, njetingap) /* Number of jets passing pfJetID and pt > 30 GeV, in pseudorapidity gap between the jets */ \
    VAR(Float_t, njetingap20) /* Number of jets passing pfJetID and pt > 20 GeV, in pseudorapidity gap between the jets */ \
    VAR(Float_t, jdphi) /* delta phi between leading jet and subleading jet */ \
    /* additional jets */ \
    VAR(Int_t, nbtag) /* pt>20 and abs(eta)<2.4 */ \
    VAR(Int_t, njets) /* pt>30 and abs(eta)<4.7 */ \
    VAR(Int_t, njetspt20) /* pt>20 and abs(eta)<4.7 */ \
    VAR(Int_t, njets_vbf) /* pt>10 and no eta cut for vbf selection */ \
    VAR(Bool_t, isVBF) /* Event is vbf if there are 2 additional jets which satisfy the VBF selection */ \
    JET_DATA(j, 1) /* leading jet sorted by pt (Fill only if corrected jet pt > 20 GeV) */ \
    JET_DATA(j, 2) /* trailing jet sorted by pt (Fill only if corrected jet pt>20 GeV) */ \
    JET_DATA(j, vbf_1) /* leading jet sorted by pt (Fill only if corrected jet pt > 20 GeV) */ \
    JET_DATA(j, vbf_2) /* trailing jet sorted by pt (Fill only if corrected jet pt>20 GeV) */ \
    JET_DATA(b, 1) /* leading b-jet sorted by pt (Fill only if corrected b-jet pt>20 GeV) */ \
    JET_DATA(b, 2) /* leading b-jet sorted by pt (Fill only if corrected b-jet pt>20 GeV) */ \
    /* Extra lepton vetos */ \
    VAR(Bool_t, extraelec_veto) /* Event is vetoed by the extra electron veto if true */ \
    VAR(Bool_t, extramuon_veto) /* Event is vetoed by the extra muon veto if true */ \
    VAR(Float_t, puweight) \
    /* hh->bbtautau part */ \
    VAR(Float_t, shapeWeight) /* genWeight * puWeight * genEventSpec */ \
    VAR(Float_t, topWeight) /* gen top pt weight for TTbar */ \
    VAR(Float_t, btagWeight) /* b tag weight */ \
    VAR(Int_t, lhe_n_partons) \
    VAR(Int_t, lhe_n_b_partons) \
    VAR(Float_t, lhe_HT) \
    VAR(Int_t, nbjets) /* pt>30 and abs(eta)<2.4 */ \
    JET_DATA(bjet_, 1) /* leading b-jet sorted by csv (Fill only if corrected b-jet pt>20 GeV) */ \
    JET_DATA(bjet_, 2) /* leading b-jet sorted by csv (Fill only if corrected b-jet pt>20 GeV) */ \
    VAR(Double_t, ht_other_jets) /* Ht of all jets in the event except the first 2 jets */\
    VAR(Float_t, m_kinfit) \
    VAR(Float_t, m_kinfit_tau_ES_up) \
    VAR(Float_t, m_kinfit_tau_ES_down) \
    VAR(Float_t, m_kinfit_jet_ES_up) \
    VAR(Float_t, m_kinfit_jet_ES_down ) \
    VAR(Int_t, kinfit_convergence) \
    VAR(Float_t, deltaR_ll) \
    VAR(UInt_t, nFatJets) \
    VAR(Int_t, hasFatJet) \
    VAR(Float_t, fatJet_pt) \
    VAR(Float_t, fatJet_eta) \
    VAR(Float_t, fatJet_phi) \
    VAR(Float_t, fatJet_energy) \
    VAR(Float_t, fatJet_m_pruned) \
    VAR(Float_t, fatJet_m_filtered) \
    VAR(Float_t, fatJet_m_trimmed) \
    VAR(Float_t, fatJet_m_softDrop) \
    VAR(Int_t, fatJet_n_subjets) \
    VAR(Float_t, fatJet_n_subjettiness_tau1) \
    VAR(Float_t, fatJet_n_subjettiness_tau2) \
    VAR(Float_t, fatJet_n_subjettiness_tau3) \
    VAR(UInt_t, genJets_nTotal) \
    VAR(UInt_t, genJets_nStored) \
    VAR(UInt_t, genJets_nStored_hadronFlavour_b) \
    VAR(UInt_t, genJets_nStored_hadronFlavour_c) \
    VAR(UInt_t, jets_nTotal_hadronFlavour_b) \
    VAR(UInt_t, jets_nTotal_hadronFlavour_c) \
    VAR(UInt_t, jets_nSelected_hadronFlavour_b) \
    VAR(UInt_t, jets_nSelected_hadronFlavour_c) \
    /**/

#define VAR(type, name) DECLARE_BRANCH_VARIABLE(type, name)
DECLARE_TREE(htt_sync, SyncEvent, SyncTuple, SYNC_DATA, "events")
#undef VAR

#define VAR(type, name) ADD_DATA_TREE_BRANCH(name)
INITIALIZE_TREE(htt_sync, SyncTuple, SYNC_DATA)
#undef VAR
#undef SYNC_DATA
#undef LEG_DATA
#undef LVAR
#undef JET_DATA
#undef JVAR

#define COND_VAL(cond, val) cond ? static_cast<float>(val) : default_value

namespace htt_sync {

    void FillSyncTuple(analysis::EventInfoBase& event, htt_sync::SyncTuple& sync, analysis::Period run_period,
                       analysis::EventInfoBase* event_tau_up = nullptr,
                       analysis::EventInfoBase* event_tau_down = nullptr,
                       analysis::EventInfoBase* event_jet_up = nullptr,
                       analysis::EventInfoBase* event_jet_down = nullptr)
    {

        static constexpr float default_value = std::numeric_limits<float>::lowest();
        static constexpr int default_int_value = std::numeric_limits<int>::lowest();

        sync().run = event->run;
        sync().lumi = event->lumi;
        sync().evt = event->evt;
        sync().sampleId = event->file_desc_id;
        sync().npv = event->npv;
        sync().npu = event->npu;

        sync().pt_1 = event->p4_1.Pt();
        sync().pt_tau_ES_up_1 = COND_VAL(event_tau_up, (*event_tau_up)->p4_1.Pt());
        sync().pt_tau_ES_down_1 = COND_VAL(event_tau_down, (*event_tau_down)->p4_1.Pt());
        sync().phi_1 = event->p4_1.Phi();
        sync().eta_1 = event->p4_1.Eta();
        sync().m_1 = event->p4_1.mass();
        sync().q_1 = event->q_1;
        sync().d0_1 = event->dxy_1;
        sync().dZ_1 = event->dz_1;
        sync().pfmt_1 = static_cast<float>(analysis::Calculate_MT(event->p4_1, event->pfMET_p4));
        sync().iso_1 =  event->iso_1;
        sync().gen_match_1 = event->gen_match_1;

        const analysis::TauIdResults tauId_1(event->tauId_flags_1);
        sync().againstElectronLooseMVA6_1 = tauId_1.Result("againstElectronLooseMVA6");
        sync().againstElectronMediumMVA6_1 = tauId_1.Result("againstElectronMediumMVA6");
        sync().againstElectronTightMVA6_1 = tauId_1.Result("againstElectronTightMVA6");
        sync().againstElectronVLooseMVA6_1 = tauId_1.Result("againstElectronVLooseMVA6");
        sync().againstElectronVTightMVA6_1 = tauId_1.Result("againstElectronVTightMVA6");
        sync().againstMuonLoose3_1 = tauId_1.Result("againstMuonLoose3");
        sync().againstMuonTight3_1 = tauId_1.Result("againstMuonTight3");
        sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = event->tauId_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
        sync().byIsolationMVArun2v1DBoldDMwLTraw_1 = event->tauId_byIsolationMVArun2v1DBoldDMwLTraw_1;
        sync().byIsolationMVArun2017v2DBoldDMwLTraw2017_1 = event->tauId_byIsolationMVArun2017v2DBoldDMwLTraw2017_1;

        sync().pt_2 = event->p4_2.Pt();
        sync().pt_tau_ES_up_2 = COND_VAL(event_tau_up, (*event_tau_up)->p4_2.Pt());
        sync().pt_tau_ES_down_2 = COND_VAL(event_tau_down, (*event_tau_down)->p4_2.Pt());
        sync().phi_2 = event->p4_2.Phi();
        sync().eta_2 = event->p4_2.Eta();
        sync().m_2 = event->p4_2.mass();
        sync().q_2 = event->q_2;
        sync().d0_2 = event->dxy_2;
        sync().dZ_2 = event->dz_2;
        sync().pfmt_2 = static_cast<float>(analysis::Calculate_MT(event->p4_2, event->pfMET_p4));
        sync().iso_2 =  event->iso_2;
        sync().gen_match_2 = event->gen_match_2;

        const analysis::TauIdResults tauId_2(event->tauId_flags_2);
        sync().againstElectronLooseMVA6_2 = tauId_2.Result("againstElectronLooseMVA6");
        sync().againstElectronMediumMVA6_2 = tauId_2.Result("againstElectronMediumMVA6");
        sync().againstElectronTightMVA6_2 = tauId_2.Result("againstElectronTightMVA6");
        sync().againstElectronVLooseMVA6_2 = tauId_2.Result("againstElectronVLooseMVA6");
        sync().againstElectronVTightMVA6_2 = tauId_2.Result("againstElectronVTightMVA6");
        sync().againstMuonLoose3_2 = tauId_2.Result("againstMuonLoose3");
        sync().againstMuonTight3_2 = tauId_2.Result("againstMuonTight3");
        sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = event->tauId_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
        sync().byIsolationMVArun2v1DBoldDMwLTraw_2 = event->tauId_byIsolationMVArun2v1DBoldDMwLTraw_2;
        sync().byIsolationMVArun2017v2DBoldDMwLTraw2017_2 = event->tauId_byIsolationMVArun2017v2DBoldDMwLTraw2017_2;

        sync().pt_tt = (event->p4_1 + event->p4_2 + event->pfMET_p4).Pt();
        sync().m_vis = (event->p4_1 + event->p4_2).M();
        sync().m_sv = event->SVfit_p4.M();
        sync().m_sv_tau_ES_up = COND_VAL(event_tau_up, (*event_tau_up)->SVfit_p4.M());
        sync().m_sv_tau_ES_down = COND_VAL(event_tau_down, (*event_tau_down)->SVfit_p4.M());
        sync().m_sv_jet_ES_up = COND_VAL(event_jet_up, (*event_jet_up)->SVfit_p4.M());
        sync().m_sv_jet_ES_down = COND_VAL(event_jet_down, (*event_jet_down)->SVfit_p4.M());
        sync().mt_sv = event->SVfit_mt;

        sync().met = event->pfMET_p4.Pt();
        sync().met_tau_ES_up = COND_VAL(event_tau_up, (*event_tau_up)->pfMET_p4.Pt());
        sync().met_tau_ES_down = COND_VAL(event_tau_down, (*event_tau_down)->pfMET_p4.Pt());
        sync().met_jet_ES_up = COND_VAL(event_jet_up, event_jet_up->GetMET().GetMomentum().Pt());
        sync().met_jet_ES_down = COND_VAL(event_jet_down, event_jet_down->GetMET().GetMomentum().Pt());

        sync().mt2 = COND_VAL(event.HasBjetPair(), event.GetMT2());
        sync().mt2_tau_ES_up = COND_VAL(event_tau_up && event_tau_up->HasBjetPair(), event_tau_up->GetMT2());
        sync().mt2_tau_ES_down = COND_VAL(event_tau_down && event_tau_down->HasBjetPair(), event_tau_down->GetMT2());
        sync().mt2_jet_ES_up = COND_VAL(event_jet_up && event_jet_up->HasBjetPair(), event_jet_up->GetMT2());
        sync().mt2_jet_ES_down = COND_VAL(event_jet_down && event_jet_down->HasBjetPair(), event_jet_down->GetMT2());

        sync().metphi = static_cast<float>(TVector2::Phi_0_2pi(event->pfMET_p4.Phi()));
        sync().metphi_tau_ES_up = COND_VAL(event_tau_up, TVector2::Phi_0_2pi((*event_tau_up)->pfMET_p4.Phi()));
        sync().metphi_tau_ES_down = COND_VAL(event_tau_down, TVector2::Phi_0_2pi((*event_tau_down)->pfMET_p4.Phi()));
        sync().metphi_jet_ES_up = COND_VAL(event_jet_up,
                                           TVector2::Phi_0_2pi(event_jet_up->GetMET().GetMomentum().Phi()));
        sync().metphi_jet_ES_down = COND_VAL(event_jet_down,
                                             TVector2::Phi_0_2pi(event_jet_down->GetMET().GetMomentum().Phi()));
        sync().pzetavis = static_cast<float>(analysis::Calculate_visiblePzeta(event->p4_1, event->p4_2));

        sync().metcov00 = static_cast<float>(event->pfMET_cov[0][0]);
        sync().metcov01 = static_cast<float>(event->pfMET_cov[0][1]);
        sync().metcov10 = static_cast<float>(event->pfMET_cov[1][0]);
        sync().metcov11 = static_cast<float>(event->pfMET_cov[1][1]);

        analysis::EventInfoBase::JetCollection jets_pt20;
        analysis::EventInfoBase::JetCollection jets_pt30;
        analysis::JetOrdering jet_ordering;
        analysis::EventInfoBase::JetCollection jets_vbf;
//        analysis::EventInfoBase::JetPair vbf_jet_pair;
        analysis::EventInfoBase::JetCollection bjets_pt;
//        analysis::EventInfoBase::JetCollection bjets_id;

        auto select_jets = [&](analysis::EventInfoBase* event_info) {
            jets_pt20.clear();
            jets_pt30.clear();
            jets_vbf.clear();
//            vbf_jet_pair = ntuple::UndefinedJetPair();
            bjets_pt.clear();
//            bjets_id.clear();
            if(!event_info) return;

            if (run_period == analysis::Period::Run2016) {
                jets_pt20 = event_info->SelectJets(20, 4.7, std::numeric_limits<double>::lowest(),
                                                   analysis::JetOrdering::Pt);
                jets_pt30 = event_info->SelectJets(30, 4.7, std::numeric_limits<double>::lowest(),
                                                   analysis::JetOrdering::Pt);
                bjets_pt = event_info->SelectJets(cuts::btag_2016::pt, cuts::btag_2016::eta, cuts::btag_2016::CSVv2M,
                                                  analysis::JetOrdering::Pt);
                bjets_id = event_info->SelectJets(cuts::btag_2016::pt, cuts::btag_2016::eta,
                                                  std::numeric_limits<double>::lowest(), analysis::JetOrdering::CSV);
                jet_ordering = analysis::JetOrdering::CSV;
            }

            if (run_period == analysis::Period::Run2017) {
                jets_pt20 = event_info->SelectJets(20, 4.7, std::numeric_limits<double>::lowest(),
                                                   analysis::JetOrdering::Pt);
                jets_pt30 = event_info->SelectJets(30, 4.7, std::numeric_limits<double>::lowest(),
                                                   analysis::JetOrdering::Pt);
                jets_vbf = event_info->SelectJets(30, 5, std::numeric_limits<double>::lowest(),
                                                  analysis::JetOrdering::Pt, event_info->GetSelectedBjetIndicesSet());
//                vbf_jet_pair = event_info->SelectVBFJetPair(jets_vbf);
                bjets_pt = event_info->SelectJets(cuts::btag_2017::pt, cuts::btag_2017::eta, cuts::btag_2017::CSVv2M,
                                                  analysis::JetOrdering::Pt);
                bjets_id = event_info->SelectJets(cuts::btag_2017::pt, cuts::btag_2017::eta,
                                                  std::numeric_limits<double>::lowest(),
                                                  analysis::JetOrdering::DeepCSV);
                jet_ordering = analysis::JetOrdering::DeepCSV;
            }

            const analysis::EventInfoBase::SelectedSignalJets selected_signal_jets = event_info->SelectSignalJets(event,run_period,jet_ordering);
        };

        select_jets(&event);

        sync().mjj = COND_VAL(jets_pt20.size() >= 2, (jets_pt20.at(0).GetMomentum()
                   + jets_pt20.at(1).GetMomentum()).M());
        sync().jdeta = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(0).GetMomentum().Eta()
                     - jets_pt20.at(1).GetMomentum().Eta());
        sync().jdphi = COND_VAL(jets_pt20.size() >= 2, TVector2::Phi_mpi_pi(jets_pt20.at(0).GetMomentum().Phi()
                     - jets_pt20.at(1).GetMomentum().Phi()));

        sync().nbtag = static_cast<int>(bjets_pt.size());
        sync().njets = static_cast<int>(jets_pt30.size());
        sync().njetspt20 = static_cast<int>(jets_pt20.size());

        sync().jpt_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
        sync().jeta_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Eta());
        sync().jrawf_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Phi());
        sync().jrawf_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0)->rawf());
        sync().jpt_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
        sync().jeta_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Eta());
        sync().jphi_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Phi());
        sync().jrawf_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1)->rawf());

        sync().njets_vbf = static_cast<int>(jets_vbf.size());
        sync().isVBF = jets_vbf.size() >= 2;
        sync().jpt_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                    jets_vbf.at(vbf_jet_pair.first).GetMomentum().Pt());
        sync().jeta_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                     jets_vbf.at(vbf_jet_pair.first).GetMomentum().Eta());
        sync().jphi_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                     jets_vbf.at(vbf_jet_pair.first).GetMomentum().Phi());
        sync().jpt_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                    jets_vbf.at(vbf_jet_pair.second).GetMomentum().Pt());
        sync().jeta_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                     jets_vbf.at(vbf_jet_pair.second).GetMomentum().Eta());
        sync().jphi_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                     jets_vbf.at(vbf_jet_pair.second).GetMomentum().Phi());

        sync().extramuon_veto = event->extramuon_veto;
        sync().extraelec_veto = event->extraelec_veto;

        sync().nbjets = static_cast<int>(bjets_id.size());
        sync().bjet_pt_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Pt());
        sync().bjet_eta_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Eta());
        sync().bjet_phi_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Phi());
        sync().bjet_rawf_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0)->rawf());
        sync().bjet_csv_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0)->csv());
        sync().bjet_deepcsv_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0)->deepcsv());
        sync().bjet_resolution_1 = COND_VAL(bjets_id.size() >= 1,
                                            bjets_id.at(0)->resolution() * bjets_id.at(0).GetMomentum().E());
        sync().bjet_pt_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Pt());
        sync().bjet_eta_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Eta());
        sync().bjet_phi_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Phi());
        sync().bjet_rawf_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1)->rawf());
        sync().bjet_csv_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1)->csv());
        sync().bjet_deepcsv_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1)->deepcsv());
        sync().bjet_resolution_2 = COND_VAL(bjets_id.size() >= 2,
                                            bjets_id.at(1)->resolution() * bjets_id.at(1).GetMomentum().E());

        sync().ht_other_jets = event.GetHT();

        sync().kinfit_convergence = COND_VAL(event->kinFit_convergence.size() > 0 && event.HasBjetPair(), event.GetKinFitResults().convergence);
        sync().m_kinfit = COND_VAL(event.GetKinFitResults().HasValidMass() && event->kinFit_convergence.size() > 0 && event.HasBjetPair(), static_cast<float>(event.GetKinFitResults().mass));
        sync().m_kinfit_tau_ES_up = COND_VAL(event_tau_up &&
                                             event_tau_up->GetKinFitResults().HasValidMass() &&
                                             event_tau_up->GetKinFitResults().convergence > 0 &&
                                             event_tau_up->HasBjetPair(), event_tau_up->GetKinFitResults().mass);
        sync().m_kinfit_tau_ES_down = COND_VAL(event_tau_down &&
                                               event_tau_down->GetKinFitResults().HasValidMass() &&
                                               event_tau_down->GetKinFitResults().convergence > 0 &&
                                               event_tau_down->HasBjetPair(), event_tau_down->GetKinFitResults().mass);
        sync().m_kinfit_jet_ES_up = COND_VAL(event_jet_up &&
                                             event_jet_up->GetKinFitResults().HasValidMass() &&
                                             event_jet_up->GetKinFitResults().convergence > 0 &&
                                             event_jet_up->HasBjetPair(), event_jet_up->GetKinFitResults().mass);
        sync().m_kinfit_jet_ES_down = COND_VAL(event_jet_down &&
                                               event_jet_down->GetKinFitResults().HasValidMass() &&
                                               event_jet_down->GetKinFitResults().convergence > 0 &&
                                               event_jet_down->HasBjetPair(), event_jet_down->GetKinFitResults().mass);



        sync().deltaR_ll = ROOT::Math::VectorUtil::DeltaR(event->p4_1, event->p4_2);

        sync().nFatJets = static_cast<unsigned>(event.GetFatJets().size());
        const auto fatJet = event.SelectFatJet(30, 0.4);
        sync().hasFatJet = bjets_id.size() >= 2 ? fatJet != nullptr : -1;
        sync().fatJet_pt = COND_VAL(fatJet, fatJet->GetMomentum().Pt());
        sync().fatJet_eta = COND_VAL(fatJet, fatJet->GetMomentum().Eta());
        sync().fatJet_phi = COND_VAL(fatJet, fatJet->GetMomentum().Phi());
        sync().fatJet_energy = COND_VAL(fatJet, fatJet->GetMomentum().E());
        sync().fatJet_m_pruned = COND_VAL(fatJet, (*fatJet)->m(ntuple::TupleFatJet::MassType::Pruned));
        sync().fatJet_m_softDrop = COND_VAL(fatJet, (*fatJet)->m(ntuple::TupleFatJet::MassType::SoftDrop));
        sync().fatJet_n_subjettiness_tau1 = COND_VAL(fatJet, (*fatJet)->n_subjettiness(1));
        sync().fatJet_n_subjettiness_tau2 = COND_VAL(fatJet, (*fatJet)->n_subjettiness(2));
        sync().fatJet_n_subjettiness_tau3 = COND_VAL(fatJet, (*fatJet)->n_subjettiness(3));
        sync().fatJet_n_subjets = fatJet ? static_cast<int>((*fatJet)->subJets().size()) : default_int_value;

        sync().topWeight = static_cast<Float_t>(event->weight_top_pt);
        sync().shapeWeight = static_cast<Float_t>(event->weight_pu * event->weight_bsm_to_sm * event->weight_dy
                           * event->weight_ttbar * event->weight_wjets * event->weight_xs);
        sync().btagWeight = static_cast<Float_t>(event->weight_btag);

        sync().lhe_n_b_partons = static_cast<int>(event->lhe_n_b_partons);
        sync().lhe_n_partons = static_cast<int>(event->lhe_n_partons);
        sync().lhe_HT = event->lhe_HT;

        sync().genJets_nTotal = event->genJets_nTotal;
        sync().genJets_nStored = static_cast<unsigned>(event->genJets_p4.size());
        sync().genJets_nStored_hadronFlavour_b = std::min<unsigned>(2, static_cast<unsigned>(
                    std::count(event->genJets_hadronFlavour.begin(), event->genJets_hadronFlavour.end(), 5)));
        sync().genJets_nStored_hadronFlavour_c = static_cast<unsigned>(
                    std::count(event->genJets_hadronFlavour.begin(), event->genJets_hadronFlavour.end(), 4));
        sync().jets_nTotal_hadronFlavour_b = event->jets_nTotal_hadronFlavour_b;
        sync().jets_nTotal_hadronFlavour_c = event->jets_nTotal_hadronFlavour_c;
        sync().jets_nSelected_hadronFlavour_b = static_cast<unsigned>(
                    std::count(event->jets_hadronFlavour.begin(), event->jets_hadronFlavour.end(), 5));
        sync().jets_nSelected_hadronFlavour_c = static_cast<unsigned>(
                    std::count(event->jets_hadronFlavour.begin(), event->jets_hadronFlavour.end(), 4));

        select_jets(event_tau_up);
        sync().jpt_tau_ES_up_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
        sync().jpt_tau_ES_up_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
        sync().bjet_pt_tau_ES_up_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Pt());
        sync().bjet_pt_tau_ES_up_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Pt());
        sync().jpt_tau_ES_up_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                              jets_vbf.at(vbf_jet_pair.first).GetMomentum().Pt());
        sync().jpt_tau_ES_up_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                              jets_vbf.at(vbf_jet_pair.second).GetMomentum().Pt());

        select_jets(event_tau_down);
        sync().jpt_tau_ES_down_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
        sync().jpt_tau_ES_down_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
        sync().bjet_pt_tau_ES_down_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Pt());
        sync().bjet_pt_tau_ES_down_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Pt());
        sync().jpt_tau_ES_down_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                                jets_vbf.at(vbf_jet_pair.first).GetMomentum().Pt());
        sync().jpt_tau_ES_down_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                                jets_vbf.at(vbf_jet_pair.second).GetMomentum().Pt());

        select_jets(event_jet_up);
        sync().jpt_jet_ES_up_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
        sync().jpt_jet_ES_up_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
        sync().bjet_pt_jet_ES_up_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Pt());
        sync().bjet_pt_jet_ES_up_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Pt());
        sync().jpt_jet_ES_up_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                              jets_vbf.at(vbf_jet_pair.first).GetMomentum().Pt());
        sync().jpt_jet_ES_up_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                              jets_vbf.at(vbf_jet_pair.second).GetMomentum().Pt());

        select_jets(event_jet_down);
        sync().jpt_jet_ES_down_1 = COND_VAL(jets_pt20.size() >= 1, jets_pt20.at(0).GetMomentum().Pt());
        sync().jpt_jet_ES_down_2 = COND_VAL(jets_pt20.size() >= 2, jets_pt20.at(1).GetMomentum().Pt());
        sync().bjet_pt_jet_ES_down_1 = COND_VAL(bjets_id.size() >= 1, bjets_id.at(0).GetMomentum().Pt());
        sync().bjet_pt_jet_ES_down_2 = COND_VAL(bjets_id.size() >= 2, bjets_id.at(1).GetMomentum().Pt());
        sync().jpt_jet_ES_down_vbf_1 = COND_VAL(vbf_jet_pair.first < jets_vbf.size(),
                                                jets_vbf.at(vbf_jet_pair.first).GetMomentum().Pt());
        sync().jpt_jet_ES_down_vbf_2 = COND_VAL(vbf_jet_pair.second < jets_vbf.size(),
                                                jets_vbf.at(vbf_jet_pair.second).GetMomentum().Pt());

        sync.Fill();
    }
}

#undef COND_VAL
