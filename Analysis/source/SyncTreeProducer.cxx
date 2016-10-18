/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/SyncTupleHTT.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/Htautau_2015.h"

struct Arguments {
    REQ_ARG(std::string, mode);
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
    REQ_ARG(std::string, output_file);
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

    SyncTreeProducer(const Arguments& _args) : args(_args)
    {
        std::istringstream ss_mode(args.mode());
        ss_mode >> syncMode;
    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% mode.\n")
                   % args.input_file() % args.output_file() % args.mode();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        EventTuple originalTuple(args.tree_name(), originalFile.get(), true, { "lhe_n_partons", "lhe_HT" });
        SyncTuple sync(args.tree_name(), outputFile.get(), false );
        const Long64_t n_entries = originalTuple.GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple.GetEntry(current_entry);
            const auto bjet_pair = EventInfoBase::SelectBjetPair(originalTuple.data(), cuts::Htautau_2015::btag::pt,
                                                                 cuts::Htautau_2015::btag::eta, JetOrdering::CSV);
            EventInfoBase event(originalTuple.data(), bjet_pair);
            if(event.GetEnergyScale() != EventEnergyScale::Central) continue;

            if(syncMode == SyncMode::HH) {
                if(event->dilepton_veto || event->extraelec_veto || event->extramuon_veto) continue;
            }

            sync().run = event->run;
            sync().lumi = event->lumi;
            sync().evt = event->evt;
            // sync().rho = ;
            sync().npv = event->npv;
            sync().npu = event->npu;

            sync().pt_1 = event->p4_1.Pt();
            sync().phi_1 = event->p4_1.Phi();
            sync().eta_1 = event->p4_1.Eta();
            sync().m_1 = event->p4_1.mass();
            sync().q_1 = event->q_1;
            sync().d0_1 = event->d0_1;
            sync().dZ_1 = event->dZ_1;
            sync().mt_1 = Calculate_MT(event->p4_1, event->mvaMET_p4);
            sync().pfmt_1 = Calculate_MT(event->p4_1, event->pfMET_p4);
            sync().puppimt_1 = Calculate_MT(event->p4_1, event->pfMET_p4);
            sync().iso_1 =  event->iso_1;
            sync().id_e_mva_nt_loose_1 = event->id_e_mva_nt_loose_1;
            sync().gen_match_1 = event->gen_match_1;
            sync().againstElectronLooseMVA6_1 = GetTauID(event->tauIDs_1, "againstElectronLooseMVA6");
            sync().againstElectronMediumMVA6_1 = GetTauID(event->tauIDs_1, "againstElectronMediumMVA6");
            sync().againstElectronTightMVA6_1 = GetTauID(event->tauIDs_1, "againstElectronTightMVA6");
            sync().againstElectronVLooseMVA6_1 = GetTauID(event->tauIDs_1, "againstElectronVLooseMVA6");
            sync().againstElectronVTightMVA6_1 = GetTauID(event->tauIDs_1, "againstElectronVTightMVA6");
            sync().againstMuonLoose3_1 = GetTauID(event->tauIDs_1, "againstMuonLoose3");
            sync().againstMuonTight3_1 = GetTauID(event->tauIDs_1, "againstMuonTight3");
            sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = GetTauID(event->tauIDs_1, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
            sync().byIsolationMVA3newDMwoLTraw_1 = GetTauID(event->tauIDs_1, "byIsolationMVA3newDMwoLTraw");
            sync().byIsolationMVA3oldDMwoLTraw_1 = GetTauID(event->tauIDs_1, "byIsolationMVA3oldDMwoLTraw");
            sync().byIsolationMVA3newDMwLTraw_1 = GetTauID(event->tauIDs_1, "byIsolationMVA3newDMwLTraw");
            sync().byIsolationMVA3oldDMwLTraw_1 = GetTauID(event->tauIDs_1, "byIsolationMVA3oldDMwLTraw");
            //sync().chargedIsoPtSum_1 = ;
            sync().decayModeFindingOldDMs_1 = GetTauID(event->tauIDs_1, "decayModeFindingOldDMs");
            // sync().neutralIsoPtSum_1 = ;
            // sync().puCorrPtSum_1 = ;
            // sync().trigweight_1 = ;
            // sync().idisoweight_1 = ;

            sync().pt_2 = event->p4_2.Pt();
            sync().phi_2 = event->p4_2.Phi();
            sync().eta_2 = event->p4_2.Eta();
            sync().m_2 = event->p4_2.mass();
            sync().q_2 = event->q_2;
            sync().d0_2 = event->d0_2;
            sync().dZ_2 = event->dZ_2;
            sync().mt_2 = Calculate_MT(event->p4_2, event->mvaMET_p4);
            sync().pfmt_2 = Calculate_MT(event->p4_2, event->pfMET_p4);
            sync().puppimt_2 = Calculate_MT(event->p4_2, event->pfMET_p4);
            sync().iso_2 =  event->iso_2;
            sync().id_e_mva_nt_loose_2 = event->id_e_mva_nt_loose_2;
            sync().gen_match_2 = event->gen_match_2;
            sync().againstElectronLooseMVA6_2 = GetTauID(event->tauIDs_2, "againstElectronLooseMVA6");
            sync().againstElectronMediumMVA6_2 = GetTauID(event->tauIDs_2, "againstElectronMediumMVA6");
            sync().againstElectronTightMVA6_2 = GetTauID(event->tauIDs_2, "againstElectronTightMVA6");
            sync().againstElectronVLooseMVA6_2 = GetTauID(event->tauIDs_2, "againstElectronVLooseMVA6");
            sync().againstElectronVTightMVA6_2 = GetTauID(event->tauIDs_2, "againstElectronVTightMVA6");
            sync().againstMuonLoose3_2 = GetTauID(event->tauIDs_2, "againstMuonLoose3");
            sync().againstMuonTight3_2 = GetTauID(event->tauIDs_2, "againstMuonTight3");
            sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = GetTauID(event->tauIDs_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
            sync().byIsolationMVA3newDMwoLTraw_2 = GetTauID(event->tauIDs_2, "byIsolationMVA3newDMwoLTraw");
            sync().byIsolationMVA3oldDMwoLTraw_2 = GetTauID(event->tauIDs_2, "byIsolationMVA3oldDMwoLTraw");
            sync().byIsolationMVA3newDMwLTraw_2 = GetTauID(event->tauIDs_2, "byIsolationMVA3newDMwLTraw");
            sync().byIsolationMVA3oldDMwLTraw_2 = GetTauID(event->tauIDs_2, "byIsolationMVA3oldDMwLTraw");
            //sync().chargedIsoPtSum_2 = ;
            sync().decayModeFindingOldDMs_2 = GetTauID(event->tauIDs_2, "decayModeFindingOldDMs");
            // sync().neutralIsoPtSum_2 = ;
            // sync().puCorrPtSum_2 = ;
            // sync().trigweight_2 = ;
            // sync().idisoweight_2 = ;

            sync().pt_tt = (event->p4_1 + event->p4_2 + event->pfMET_p4).Pt();
            sync().mt_tot = Calculate_TotalMT(event->p4_1, event->p4_2, event->mvaMET_p4);
            sync().m_vis = (event->p4_1 + event->p4_2).M();
            sync().m_sv = event->SVfit_p4.M();
            //sync().mt_sv = ;

            sync().met = event->pfMET_p4.Pt();
            sync().metphi = TVector2::Phi_0_2pi(event->pfMET_p4.Phi());
            sync().puppimet = event->puppiMET_p4.Pt();
            sync().puppimetphi = event->puppiMET_p4.Phi();
            sync().mvamet = event->mvaMET_p4.Pt();
            sync().mvametphi = event->mvaMET_p4.Phi();
            sync().pzetavis = Calculate_visiblePzeta(event->p4_1, event->p4_2);
            sync().pzetamiss = Calculate_Pzeta(event->p4_1, event->p4_2, event->mvaMET_p4);
            sync().mvacov00 = event->mvaMET_cov[0][0];
            sync().mvacov01 = event->mvaMET_cov[0][1];
            sync().mvacov10 = event->mvaMET_cov[1][0];
            sync().mvacov11 = event->mvaMET_cov[1][1];
            sync().metcov00 = event->pfMET_cov[0][0];
            sync().metcov01 = event->pfMET_cov[0][1];
            sync().metcov10 = event->pfMET_cov[1][0];
            sync().metcov11 = event->pfMET_cov[1][1];

            const auto jets_pt20 = event.SelectJets(20, 4.7, std::numeric_limits<double>::lowest(), JetOrdering::Pt);
            const auto jets_pt30 = event.SelectJets(30, 4.7, std::numeric_limits<double>::lowest(), JetOrdering::Pt);
            const auto bjets_pt = event.SelectJets(20, 2.4, std::numeric_limits<double>::lowest(), JetOrdering::Pt);
            const auto bjets_csv = event.SelectJets(20, 2.4, std::numeric_limits<double>::lowest(), JetOrdering::CSV);

            if(jets_pt20.size() >= 2) {
                sync().mjj = (jets_pt20.at(0).GetMomentum() + jets_pt20.at(1).GetMomentum()).M();
                sync().jdeta = jets_pt20.at(0).GetMomentum().Eta() - jets_pt20.at(1).GetMomentum().Eta();
                //sync().njetingap = ;
                //sync().njetingap20 = ;
                sync().jdphi = TVector2::Phi_mpi_pi(jets_pt20.at(0).GetMomentum().Phi()
                                                    - jets_pt20.at(1).GetMomentum().Phi());
            } else {
                sync().mjj = default_value;
                sync().jdeta = default_value;
                //sync().njetingap = default_value;
                //sync().njetingap20 = default_value;
                sync().jdphi = default_value;
            }

            sync().nbtag = bjets_pt.size();
            sync().njets = jets_pt30.size();
            sync().njetspt20 = jets_pt20.size();
            if(jets_pt20.size() >= 1) {
                sync().jpt_1 = jets_pt20.at(0).GetMomentum().Pt();
                sync().jeta_1 = jets_pt20.at(0).GetMomentum().Eta();
                sync().jphi_1 = jets_pt20.at(0).GetMomentum().Phi();
                sync().jrawf_1 = jets_pt20.at(0)->rawf();
                sync().jmva_1 = jets_pt20.at(0)->mva();
            } else {
                sync().jpt_1 = default_value;
                sync().jeta_1 = default_value;
                sync().jphi_1 = default_value;
                sync().jrawf_1 = default_value;
                sync().jmva_1 = default_value;
            }
            if(jets_pt20.size() >= 2) {
                sync().jpt_2 = jets_pt20.at(1).GetMomentum().Pt();
                sync().jeta_2 = jets_pt20.at(1).GetMomentum().Eta();
                sync().jphi_2 = jets_pt20.at(1).GetMomentum().Phi();
                sync().jrawf_2 = jets_pt20.at(1)->rawf();
                sync().jmva_2 = jets_pt20.at(1)->mva();
            } else {
                sync().jpt_2 = default_value;
                sync().jeta_2 = default_value;
                sync().jphi_2 = default_value;
                sync().jrawf_2 = default_value;
                sync().jmva_2 = default_value;
            }
            if(bjets_pt.size() >= 1) {
                sync().bpt_1 = bjets_pt.at(0).GetMomentum().Pt();
                sync().beta_1 = bjets_pt.at(0).GetMomentum().Eta();
                sync().bphi_1 = bjets_pt.at(0).GetMomentum().Phi();
                sync().brawf_1 = bjets_pt.at(0)->rawf();
                sync().bmva_1 = bjets_pt.at(0)->mva();
                sync().bcsv_1 = bjets_pt.at(0)->csv();
            } else {
                sync().bpt_1 = default_value;
                sync().beta_1 = default_value;
                sync().bphi_1 = default_value;
                sync().brawf_1 = default_value;
                sync().bmva_1 = default_value;
                sync().bcsv_1 = default_value;
            }
            if(bjets_pt.size() >= 2) {
                sync().bpt_2 = bjets_pt.at(1).GetMomentum().Pt();
                sync().beta_2 = bjets_pt.at(1).GetMomentum().Eta();
                sync().bphi_2 = bjets_pt.at(1).GetMomentum().Phi();
                sync().brawf_2 = bjets_pt.at(1)->rawf();
                sync().bmva_2 = bjets_pt.at(1)->mva();
                sync().bcsv_2 = bjets_pt.at(1)->csv();
            } else {
                sync().bpt_2 = default_value;
                sync().beta_2 = default_value;
                sync().bphi_2 = default_value;
                sync().brawf_2 = default_value;
                sync().bmva_2 = default_value;
                sync().bcsv_2 = default_value;
            }

            sync().dilepton_veto = event->dilepton_veto;
            sync().extramuon_veto = event->extramuon_veto;
            sync().extraelec_veto = event->extraelec_veto;
//            sync().puweight = ;

            // hh->bbtautau part
            sync().nbjets = bjets_csv.size();
            if(bjets_csv.size() >= 1) {
                sync().bjet_pt_1 = bjets_csv.at(0).GetMomentum().Pt();
                sync().bjet_eta_1 = bjets_csv.at(0).GetMomentum().Eta();
                sync().bjet_phi_1 = bjets_csv.at(0).GetMomentum().Phi();
                sync().bjet_rawf_1 = bjets_csv.at(0)->rawf();
                sync().bjet_mva_1 = bjets_csv.at(0)->mva();
                sync().bjet_csv_1 = bjets_csv.at(0)->csv();
            } else {
                sync().bjet_pt_1 = default_value;
                sync().bjet_eta_1 = default_value;
                sync().bjet_phi_1 = default_value;
                sync().bjet_rawf_1 = default_value;
                sync().bjet_mva_1 = default_value;
                sync().bjet_csv_1 = default_value;
            }
            if(bjets_csv.size() >= 2) {
                sync().bjet_pt_2 = bjets_csv.at(1).GetMomentum().Pt();
                sync().bjet_eta_2 = bjets_csv.at(1).GetMomentum().Eta();
                sync().bjet_phi_2 = bjets_csv.at(1).GetMomentum().Phi();
                sync().bjet_rawf_2 = bjets_csv.at(1)->rawf();
                sync().bjet_mva_2 = bjets_csv.at(1)->mva();
                sync().bjet_csv_2 = bjets_csv.at(1)->csv();
            } else {
                sync().bjet_pt_2 = default_value;
                sync().bjet_eta_2 = default_value;
                sync().bjet_phi_2 = default_value;
                sync().bjet_rawf_2 = default_value;
                sync().bjet_mva_2 = default_value;
                sync().bjet_csv_2 = default_value;
            }

//            if(bjets_csv.size() >= 2)
//                sync().kinfit_convergence = event.GetKinFitResults().convergence;
//            else
//                sync().kinfit_convergence = default_int_value;

//            if(bjets_csv.size() >= 2 && event.GetKinFitResults().HasValidMass())
//                sync().m_kinfit = event.GetKinFitResults().mass;
//            else
//                sync().m_kinfit = default_value;

            sync().deltaR_ll = ROOT::Math::VectorUtil::DeltaR(event->p4_1, event->p4_2);

            sync.Fill();
        }
        sync.Write();
    }

private:
    float GetTauID(const std::map<std::string, float>& id_map, const std::string& id_name) const
    {
        return id_map.count(id_name) ? id_map.at(id_name) : default_value;
    }

private:
    Arguments args;
    SyncMode syncMode;
};

} // namespace analysis

PROGRAM_MAIN(analysis::SyncTreeProducer, Arguments)
