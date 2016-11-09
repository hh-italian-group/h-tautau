/*! Produce synchronization tree.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "AnalysisTools/Run/include/program_main.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/SyncTupleHTT.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/McCorrections/include/EventWeights.h"

struct Arguments {
    REQ_ARG(std::string, mode);
    REQ_ARG(std::string, input_file);
    REQ_ARG(std::string, tree_name);
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
    }

    void Run()
    {
        std::cout << boost::format("Processing input file '%1%' into output file '%2%' using %3% mode.\n")
                   % args.input_file() % args.output_file() % args.mode();

        auto originalFile = root_ext::OpenRootFile(args.input_file());
        auto outputFile = root_ext::CreateRootFile(args.output_file());
        EventTuple originalTuple(args.tree_name(), originalFile.get(), true);
        SyncTuple sync(args.tree_name(), outputFile.get(), false );
        const Long64_t n_entries = originalTuple.GetEntries();
        for(Long64_t current_entry = 0; current_entry < n_entries; ++current_entry) {
            originalTuple.GetEntry(current_entry);
            const auto bjet_pair = EventInfoBase::SelectBjetPair(originalTuple.data(), cuts::btag_2016::pt,
                                                                 cuts::btag_2016::eta, JetOrdering::CSV);
            EventInfoBase event(originalTuple.data(), bjet_pair);
            if(event.GetEnergyScale() != EventEnergyScale::Central || !event->trigger_match) continue;

            if(syncMode == SyncMode::HH) {
                if(/*event->dilepton_veto ||*/ event->extraelec_veto || event->extramuon_veto) continue;
            }

            const auto GetTauID = [&](size_t leg_id, const std::string& id_name) -> float {
                const auto& id_keys = leg_id == 1 ? event->tauId_keys_1 : event->tauId_keys_2;
                const auto& id_values = leg_id == 1 ? event->tauId_values_1 : event->tauId_values_2;
                const uint32_t key = tools::hash(id_name);
                const auto iter = std::find(id_keys.begin(), id_keys.end(), key);
                if(iter == id_keys.end()) return default_value;
                const auto index = std::distance(id_keys.begin(), iter);
                return id_values.at(index);
            };


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
            sync().d0_1 = event->dxy_1;
            sync().dZ_1 = event->dz_1;
//            sync().mt_1 = Calculate_MT(event->p4_1, event->mvaMET_p4);
            sync().pfmt_1 = Calculate_MT(event->p4_1, event->pfMET_p4);
//            sync().puppimt_1 = Calculate_MT(event->p4_1, event->pfMET_p4);
            sync().iso_1 =  event->iso_1;
//            sync().id_e_mva_nt_loose_1 = event->id_e_mva_nt_loose_1;
            sync().gen_match_1 = event->gen_match_1;
            sync().againstElectronLooseMVA6_1 = GetTauID(1, "againstElectronLooseMVA6");
            sync().againstElectronMediumMVA6_1 = GetTauID(1, "againstElectronMediumMVA6");
            sync().againstElectronTightMVA6_1 = GetTauID(1, "againstElectronTightMVA6");
            sync().againstElectronVLooseMVA6_1 = GetTauID(1, "againstElectronVLooseMVA6");
            sync().againstElectronVTightMVA6_1 = GetTauID(1, "againstElectronVTightMVA6");
            sync().againstMuonLoose3_1 = GetTauID(1, "againstMuonLoose3");
            sync().againstMuonTight3_1 = GetTauID(1, "againstMuonTight3");
            sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_1 = GetTauID(1, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
            sync().byIsolationMVA3newDMwoLTraw_1 = GetTauID(1, "byIsolationMVA3newDMwoLTraw");
            sync().byIsolationMVA3oldDMwoLTraw_1 = GetTauID(1, "byIsolationMVA3oldDMwoLTraw");
            sync().byIsolationMVA3newDMwLTraw_1 = GetTauID(1, "byIsolationMVA3newDMwLTraw");
            sync().byIsolationMVA3oldDMwLTraw_1 = GetTauID(1, "byIsolationMVA3oldDMwLTraw");
            //sync().chargedIsoPtSum_1 = ;
            sync().decayModeFindingOldDMs_1 = GetTauID(1, "decayModeFindingOldDMs");
            // sync().neutralIsoPtSum_1 = ;
            // sync().puCorrPtSum_1 = ;
            // sync().trigweight_1 = ;
            // sync().idisoweight_1 = ;

            sync().pt_2 = event->p4_2.Pt();
            sync().phi_2 = event->p4_2.Phi();
            sync().eta_2 = event->p4_2.Eta();
            sync().m_2 = event->p4_2.mass();
            sync().q_2 = event->q_2;
            sync().d0_2 = event->dxy_2;
            sync().dZ_2 = event->dz_2;
//            sync().mt_2 = Calculate_MT(event->p4_2, event->mvaMET_p4);
            sync().pfmt_2 = Calculate_MT(event->p4_2, event->pfMET_p4);
//            sync().puppimt_2 = Calculate_MT(event->p4_2, event->pfMET_p4);
            sync().iso_2 =  event->iso_2;
//            sync().id_e_mva_nt_loose_2 = event->id_e_mva_nt_loose_2;
            sync().gen_match_2 = event->gen_match_2;
            sync().againstElectronLooseMVA6_2 = GetTauID(2, "againstElectronLooseMVA6");
            sync().againstElectronMediumMVA6_2 = GetTauID(2, "againstElectronMediumMVA6");
            sync().againstElectronTightMVA6_2 = GetTauID(2, "againstElectronTightMVA6");
            sync().againstElectronVLooseMVA6_2 = GetTauID(2, "againstElectronVLooseMVA6");
            sync().againstElectronVTightMVA6_2 = GetTauID(2, "againstElectronVTightMVA6");
            sync().againstMuonLoose3_2 = GetTauID(2, "againstMuonLoose3");
            sync().againstMuonTight3_2 = GetTauID(2, "againstMuonTight3");
            sync().byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = GetTauID(2, "byCombinedIsolationDeltaBetaCorrRaw3Hits");
            sync().byIsolationMVA3newDMwoLTraw_2 = GetTauID(2, "byIsolationMVA3newDMwoLTraw");
            sync().byIsolationMVA3oldDMwoLTraw_2 = GetTauID(2, "byIsolationMVA3oldDMwoLTraw");
            sync().byIsolationMVA3newDMwLTraw_2 = GetTauID(2, "byIsolationMVA3newDMwLTraw");
            sync().byIsolationMVA3oldDMwLTraw_2 = GetTauID(2, "byIsolationMVA3oldDMwLTraw");
            //sync().chargedIsoPtSum_2 = ;
            sync().decayModeFindingOldDMs_2 = GetTauID(2, "decayModeFindingOldDMs");
            // sync().neutralIsoPtSum_2 = ;
            // sync().puCorrPtSum_2 = ;
            // sync().trigweight_2 = ;
            // sync().idisoweight_2 = ;

            sync().pt_tt = (event->p4_1 + event->p4_2 + event->pfMET_p4).Pt();
//            sync().mt_tot = Calculate_TotalMT(event->p4_1, event->p4_2, event->mvaMET_p4);
            sync().m_vis = (event->p4_1 + event->p4_2).M();
            sync().m_sv = event->SVfit_p4.M();
            sync().mt_sv = event->SVfit_mt;

            sync().met = event->pfMET_p4.Pt();
            if(syncMode == SyncMode::HH) sync().metphi = TVector2::Phi_0_2pi(event->pfMET_p4.Phi());
            else sync().metphi = TVector2::Phi_mpi_pi(event->pfMET_p4.Phi());
//            sync().puppimet = event->puppiMET_p4.Pt();
//            sync().puppimetphi = event->puppiMET_p4.Phi();
//            sync().mvamet = event->mvaMET_p4.Pt();
//            sync().mvametphi = event->mvaMET_p4.Phi();
            sync().pzetavis = Calculate_visiblePzeta(event->p4_1, event->p4_2);
//            sync().pzetamiss = Calculate_Pzeta(event->p4_1, event->p4_2, event->mvaMET_p4);
//            sync().mvacov00 = event->mvaMET_cov[0][0];
//            sync().mvacov01 = event->mvaMET_cov[0][1];
//            sync().mvacov10 = event->mvaMET_cov[1][0];
//            sync().mvacov11 = event->mvaMET_cov[1][1];
            sync().metcov00 = event->pfMET_cov[0][0];
            sync().metcov01 = event->pfMET_cov[0][1];
            sync().metcov10 = event->pfMET_cov[1][0];
            sync().metcov11 = event->pfMET_cov[1][1];

            const auto jets_pt20 = event.SelectJets(20, 4.7, std::numeric_limits<double>::lowest(), JetOrdering::Pt);
            const auto jets_pt30 = event.SelectJets(30, 4.7, std::numeric_limits<double>::lowest(), JetOrdering::Pt);
            const auto bjets_pt = event.SelectJets(cuts::btag_2016::pt, cuts::btag_2016::eta, cuts::btag_2016::CSVv2M, JetOrdering::Pt);
            const auto bjets_csv = event.SelectJets(cuts::btag_2016::pt,cuts::btag_2016::eta,std::numeric_limits<double>::lowest(), JetOrdering::CSV);

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

            if(syncMode == SyncMode::HH){

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


                if(event->kinFit_convergence.size() > 0) {
                    if(bjets_csv.size() >= 2)
                        sync().kinfit_convergence = event.GetKinFitResults().convergence;
                    else
                        sync().kinfit_convergence = default_int_value;

                    if(bjets_csv.size() >= 2 && event.GetKinFitResults().HasValidMass())
                        sync().m_kinfit = event.GetKinFitResults().mass;
                    else
                        sync().m_kinfit = default_value;
                }

                sync().deltaR_ll = ROOT::Math::VectorUtil::DeltaR(event->p4_1, event->p4_2);

                sync().nFatJets = event.GetFatJets().size();
                const FatJetCandidate* fatJetPtr = event.SelectFatJet(30, 0.4);
                sync().hasFatJet = bjets_csv.size() >= 2 ? fatJetPtr != nullptr : -1;
                if(fatJetPtr) {
                    const FatJetCandidate& fatJet = *fatJetPtr;
                    sync().fatJet_pt = fatJet.GetMomentum().Pt();
                    sync().fatJet_eta = fatJet.GetMomentum().Eta();
                    sync().fatJet_phi = fatJet.GetMomentum().Phi();
                    sync().fatJet_energy = fatJet.GetMomentum().E();
                    sync().fatJet_m_pruned = fatJet->m(ntuple::TupleFatJet::MassType::Pruned);
                    //                sync().fatJet_m_filtered = fatJet->m(ntuple::TupleFatJet::MassType::Filtered);
                    //                sync().fatJet_m_trimmed = fatJet->m(ntuple::TupleFatJet::MassType::Trimmed);
                    sync().fatJet_m_softDrop = fatJet->m(ntuple::TupleFatJet::MassType::SoftDrop);
                    sync().fatJet_n_subjets = fatJet->subJets().size();
                    sync().fatJet_n_subjettiness_tau1 = fatJet->n_subjettiness(1);
                    sync().fatJet_n_subjettiness_tau2 = fatJet->n_subjettiness(2);
                    sync().fatJet_n_subjettiness_tau3 = fatJet->n_subjettiness(3);
                } else {
                    sync().fatJet_pt = default_value;
                    sync().fatJet_eta = default_value;
                    sync().fatJet_phi = default_value;
                    sync().fatJet_energy = default_value;
                    sync().fatJet_m_pruned = default_value;
                    sync().fatJet_m_filtered = default_value;
                    sync().fatJet_m_trimmed = default_value;
                    sync().fatJet_m_softDrop = default_value;
                    sync().fatJet_n_subjets = default_int_value;
                    sync().fatJet_n_subjettiness_tau1 = default_value;
                    sync().fatJet_n_subjettiness_tau2 = default_value;
                    sync().fatJet_n_subjettiness_tau3 = default_value;
                }

                double topWeight = 1;
                if(args.sample_type() == "ttbar") {
                    for(size_t n = 0; n < event->genParticles_pdg.size(); ++n) {
                        if(std::abs(event->genParticles_pdg.at(n)) != 6) continue;
                        const double pt = event->genParticles_p4.at(n).pt();
                        topWeight *= std::sqrt(std::exp(0.156 - 0.00137 * pt));
                    }
                }
                sync().topWeight = topWeight;
                sync().shapeWeight = eventWeights.GetPileUpWeight(*event) * event->genEventWeight;
                sync().lhe_n_b_partons = event->lhe_n_b_partons;
                sync().lhe_n_partons = event->lhe_n_partons;
                sync().lhe_HT = event->lhe_HT;
            }

            sync.Fill();
        }

        sync.Write();
    }

private:
    Arguments args;
    SyncMode syncMode;
    mc_corrections::EventWeights eventWeights;
};

} // namespace analysis

PROGRAM_MAIN(analysis::SyncTreeProducer, Arguments)
