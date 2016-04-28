/*! Implementation of a SyncTree producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Production/interface/SyncTreeProducer_tautau.h"

//
// constructors and destructor
//
SyncTreeProducer_tautau::SyncTreeProducer_tautau(const edm::ParameterSet& iConfig):
  BaseEDAnalyzer(iConfig),
  sync_tree(BaseEDAnalyzer::GetSyncTree(&edm::Service<TFileService>()->file())),
  syncTree(*sync_tree),
  anaData("tautau_cuts.root") {}


SyncTreeProducer_tautau::~SyncTreeProducer_tautau()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

analysis::Channel SyncTreeProducer_tautau::ChannelId() const {return analysis::Channel::TauTau;}

//
// member functions
//

// ------------ method called for each event  ------------
void
SyncTreeProducer_tautau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::TauTau;

  const auto Key = analysis::stringToDataSourceTypeMap.at(BaseEDAnalyzer::GetSampleType());
  const auto& hltPaths = TauTau::trigger::hltPathMaps.at(Key);

  cuts::Cutter cut(&GetAnaData().Selection("events"));

  //Get collection
  Initialize(iEvent);

  try{

      cut(true,"events");

      cut(HaveTriggerFired(iEvent, hltPaths),"trigger");

      auto lheInfoPair    = BaseEDAnalyzer::computeHtValue();
      selection.HT        = lheInfoPair.first;
      selection.weightevt = (BaseEDAnalyzer::GetGenEventInfo())->weight();

      const auto PV = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0); //Deferenzio per passare da edm::Handle a edm::View. Quest'ultimo permette
                                            //di gestire una qualsiasi collezione del tipo passatogli tramite Tamplate.
                                            //Es. edm::View<int> gestisce int semplici, vector<int>, set<int> etc.
      cut(PV.isNonnull(),"vertex");

      selection.npv = (*(BaseEDAnalyzer::GetVertexCollection())).size();

      if(PUInfo.isValid()){
              for(std::vector<PileupSummaryInfo>::const_iterator PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI){
                  int BX = PVI->getBunchCrossing();
                  if(BX == 0) selection.numtruepileupinteractions = PVI->getTrueNumInteractions();
            }
          }

//      const auto selectedElectrons = CollectSignalElectrons();

//      cut(selectedElectrons.size(),"electrons");

      const auto selectedTaus = CollectSignalTaus();

      cut(selectedTaus.size(),"taus");

//      auto higgses = BaseEDAnalyzer::FindCompatibleObjects(selectedElectrons,selectedTaus,DeltaR_betweenSignalObjects,analysis::CandidateV2::Type::Higgs,
//                                         "H_mu_tau");

//      cut(higgses.size(),"ele_tau_pair");

//      bool isCrossTrigger = false;
//      auto triggeredHiggses = ApplyTriggerMatch(iEvent, higgses,hltPaths,false,isCrossTrigger);

//      cut(triggeredHiggses.size(),"triggerMatch");

//      selection.higgs = SelectHiggs(triggeredHiggses);
//      auto higgs = selection.higgs;

//      //Third-Lepton Veto
//      const auto muonVetoCollection     = CollectVetoMuons();
//      const auto electronVetoCollection = CollectVetoElectrons();

//      selection.muonVeto = muonVetoCollection.size() ? true : false;

//      selection.electronVeto = false;
//      for (auto &electron: electronVetoCollection){
//           if (!(selection.GetLeg(1)->GetMomentum().DeltaR(electron->GetMomentum()) < 0.05)) selection.electronVeto = true;
//      }

//      analysis::MissingETPtr pfMET(new analysis::MissingET((*pfMETs)[0],*(BaseEDAnalyzer::GetMETCovMatrix())));
//      selection.pfMET = pfMET;

//      selection.jets   = CollectJets();
//      selection.bjets  = CollectBJets();

//      selection.svfitResult = BaseEDAnalyzer::SVFit<pat::Electron>(higgs,pfMET);
//          std::cout<< "\n\t CHOOSEN SVFit -->  " << selection.svfitResult.mass <<std::endl;

      FillSyncTree(iEvent);

  }catch(cuts::cut_failed&){}

  GetAnaData().Selection("events").fill_selection();

}

// ------------ objects selection  ------------

void SyncTreeProducer_tautau::SelectSignalTau(const analysis::CandidateV2Ptr& tau, analysis::SelectionManager& selectionManager, cuts::Cutter& cut){
    using namespace cuts::Htautau_2015::ETau;
    using namespace cuts::Htautau_2015::ETau::tauID;
    const pat::Tau& object = tau->GetNtupleObject<pat::Tau>();

    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(object.leadChargedHadrCand().get());

    double const taudz = packedLeadTauCand->dz();
    int const decayModeFinding    = object.tauID("decayModeFinding");

    cut(true, ">0 tau cand");
    cut(X(pt()) > pt, "pt");
    cut(std::fabs( X(eta()) ) < eta, "eta");
    cut(Y(decayModeFinding) > tauID::decayModeFinding, "oldDecayMode");
    cut(std::fabs( Y(taudz) ) < dz ,"dz");
    cut(std::abs(X(charge()))==1,"charge");
}

void SyncTreeProducer_tautau::SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_2015;
        using namespace cuts::Htautau_2015::jetID;
        const pat::Jet& object = jet->GetNtupleObject<pat::Jet>();

        cut(true, ">0 jet cand");
        cut(X(pt()) > pt_loose, "pt_loose");
        cut(std::abs( X(eta()) ) < eta, "eta");
        const bool jetMVAID = passPFLooseId(object);
        cut(Y(jetMVAID),"jet_id");
        const double deltaR_leg1 = jet->GetMomentum().DeltaR(selection.GetLeg(1)->GetMomentum());
        cut(Y(deltaR_leg1) > deltaR_signalObjects,"deltaR_lep1");
        const double deltaR_leg2 = jet->GetMomentum().DeltaR(selection.GetLeg(2)->GetMomentum());
        cut(Y(deltaR_leg2) > deltaR_signalObjects,"deltaR_lep2");
    }

void SyncTreeProducer_tautau::SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
{
    using namespace cuts::Htautau_2015;
    using namespace cuts::Htautau_2015::jetID;
    const pat::Jet& object = jet->GetNtupleObject<pat::Jet>();

    SelectJets(jet,selectionManager, cut);
    cut(std::abs( X(eta()) ) < btag::eta, "eta");
    const double csvValue = object.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    cut(Y(csvValue) > btag::CSV, "btaggin");

}


// ------------ method called once each job just before starting event loop  ------------
void
SyncTreeProducer_tautau::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SyncTreeProducer_tautau::endJob()
{
    syncTree.Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void
SyncTreeProducer_mutau::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
SyncTreeProducer_mutau::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SyncTreeProducer_mutau::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SyncTreeProducer_mutau::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//void
//SyncTreeProducer_etau::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//  //The following says we do not know what parameters are allowed so do no validation
//  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);
//}


analysis::CandidateV2Ptr SyncTreeProducer_tautau::SelectHiggs(analysis::CandidateV2PtrVector& higgses)
{
    if(!higgses.size())
        throw std::runtime_error("no available higgs candidate to select");
    if(higgses.size()==1) return higgses.front();

    const auto higgsSelector = [&] (const analysis::CandidateV2Ptr& first, const analysis::CandidateV2Ptr& second) -> bool
    {
        const pat::Electron& first_electron  = first->GetDaughter(analysis::CandidateV2::Type::Electron)->GetNtupleObject<pat::Electron>();
        const pat::Tau&  first_tau           = first->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();
        const pat::Electron& second_electron = second->GetDaughter(analysis::CandidateV2::Type::Electron)->GetNtupleObject<pat::Electron>();
        const pat::Tau&  second_tau          = second->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();

        double iso_ele1 = (first_electron.pfIsolationVariables().sumChargedHadronPt
                           + std::max(first_electron.pfIsolationVariables().sumNeutralHadronEt
                           + first_electron.pfIsolationVariables().sumPhotonEt -
                           0.5 * first_electron.pfIsolationVariables().sumPUPt, 0.0)) / first_electron.pt();
        double iso_ele2 = (second_electron.pfIsolationVariables().sumChargedHadronPt
                           + std::max(second_electron.pfIsolationVariables().sumNeutralHadronEt
                           + second_electron.pfIsolationVariables().sumPhotonEt -
                           0.5 * second_electron.pfIsolationVariables().sumPUPt, 0.0)) / second_electron.pt();

        double iso_tau1 = first_tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        double iso_tau2 = second_tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");

        std::cout << "LAMBDA ------------------------  \n IsoMu1 = "<<iso_ele1<<"  IsoMu2 = "<<iso_ele2
                  << "  IsoTau1 = "<<iso_tau1<<"  IsoTau2 = "<<iso_tau2<<std::endl;
        std::cout << " PtMu1 = "<<first_electron.pt()<<"  PtMu2 = "<<second_electron.pt()
                  << "  PtTau1 = "<<first_tau.pt()<<"  PtTau2 = "<<second_tau.pt()<<std::endl;

        bool electron = (iso_ele1 < iso_ele2) ||
                        ((iso_ele1 == iso_ele2) ? first_electron.pt() > second_electron.pt() : false);
        bool tau      = (iso_tau1 > iso_tau2) ||
                        ((iso_tau1 == iso_tau2) ? first_tau.pt() > second_tau.pt() : false);

//        if ( &first_muon==&second_muon ){
//            if ( &first_tau==&second_tau ) return true;
//            return tau;
//        }

//        if (!muon) return tau;
//        return muon;
        if ( &first_electron!=&second_electron ) return electron;
        if ( &first_tau!=&second_tau ) return tau;

        throw analysis::exception("not found a good criteria for best tau pair");
    };

    std::sort(higgses.begin(), higgses.end(), higgsSelector) ;
    return higgses.front();
}

void SyncTreeProducer_tautau::FillSyncTree(const edm::Event& iEvent)
    {
//         using namespace analysis;
//         static const float default_value = Run2::DefaultFloatFillValueForSyncTree();
//
//
//         const VertexV2Ptr primaryVertex(new VertexV2((*(BaseEDAnalyzer::GetVertexCollection())).front()));
//
//         //BaseEDAnalyzer::FillSyncTree(iEvent, selection, syncTree);
//
//         //const auto primaryVertex = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);
//
//         // Event
//         std::cout<<"~~~~~~~~~~~~~EVENT Info~~~~~~~~~"<<std::endl;
//         std::cout<<"Run = "<<iEvent.id().run()<<"  Lumi = "<<iEvent.id().luminosityBlock()
//                 <<" Event = "<<iEvent.id().event()<<std::endl;
//         syncTree.run()  = iEvent.id().run();
//         syncTree.lumi() = iEvent.id().luminosityBlock();
//         syncTree.evt()  = iEvent.id().event();
//         syncTree.channelID() = static_cast<int>(SyncTreeProducer_tautau::ChannelId());
//
// //         if(computeHT_){
// //             syncTree.HT()   = selection.HT;
// //             if(selection.HT<100) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::lt100);
// //             if(100<=selection.HT && selection.HT<200) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f100to200);
// //             if(200<=selection.HT && selection.HT<400) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f200to400);
// //             if(400<=selection.HT && selection.HT<600) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::f400to600);
// //             if(selection.HT>=600) syncTree.HTBin() = static_cast<int>(Run2::HTbinning::gt600);
// // //        syncTree->eventType() = static_cast<int>(selection.eventType);
// // //        syncTree->eventEnergyScale() = static_cast<int>(eventEnergyScale);
// //         }
// //         else {
// //             syncTree.HT() = default_value;
// //             syncTree.HTBin() = -1;
// //         }
//
//         // if (sampleType == "Spring15MC" || sampleType == "Fall15MC") syncTree.weightevt() = selection.weightevt;
// //         else syncTree.weightevt() = default_value;
//
//         syncTree.weightevt() = default_value;
//         syncTree.npv() = selection.npv;
//         syncTree.npu() = selection.numtruepileupinteractions;
//
//
//         // HTT candidate
//         syncTree.m_vis() = selection.higgs->GetMomentum().M();
//         syncTree.pt_tt()  = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()).Pt();
// //        syncTree.m_sv() = selection.svfitResult.has_valid_mass
// //                ? selection.svfitResult.mass : Run2::DefaultFillValueForSyncTree();
// //        syncTree.pt_sv() = selection.svfitResult.has_valid_momentum
// //                ? selection.svfitResult.momentum.Pt() : Run2::DefaultFillValueForSyncTree();
// //        syncTree.eta_sv() = selection.svfitResult.has_valid_momentum
// //                ? selection.svfitResult.momentum.Eta() : Run2::DefaultFillValueForSyncTree();
// //        syncTree.phi_sv() = selection.svfitResult.has_valid_momentum
// //                ? selection.svfitResult.momentum.Phi() : Run2::DefaultFillValueForSyncTree();
//
//
//         syncTree.met()      = selection.pfMET->Pt();
//         syncTree.metphi()   = selection.pfMET->Phi();
//         syncTree.isPFMET()  = selection.pfMET->isPFMET();
//         syncTree.metcov00() = selection.pfMET->GetCovVector().at(0);
//         syncTree.metcov01() = selection.pfMET->GetCovVector().at(1);
//         syncTree.metcov10() = selection.pfMET->GetCovVector().at(2);
//         syncTree.metcov11() = selection.pfMET->GetCovVector().at(3);
//
//         // Leg 1, lepton
//         const pat::Electron& patElectron = selection.GetLeg(1)->GetNtupleObject<pat::Electron>();
//
//         syncTree.pt_1()     = selection.GetLeg(1)->GetMomentum().Pt();
//         syncTree.phi_1()    = selection.GetLeg(1)->GetMomentum().Phi();
//         syncTree.eta_1()    = selection.GetLeg(1)->GetMomentum().Eta();
//         syncTree.m_1()      = selection.GetLeg(1)->GetMomentum().M();
//         syncTree.q_1()      = selection.GetLeg(1)->GetCharge();
//         syncTree.pfmt_1()     = Calculate_MT(selection.GetLeg(1)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
//         syncTree.d0_1()     = Calculate_dxy(selection.GetLeg(1)->GetVertexPosition(), primaryVertex->GetPosition(),
//                                              selection.GetLeg(1)->GetMomentum());
//         syncTree.dZ_1()     = selection.GetLeg(1)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
//         syncTree.iso_1()    = (patElectron.pfIsolationVariables().sumChargedHadronPt
//                                + std::max(patElectron.pfIsolationVariables().sumNeutralHadronEt
//                                + patElectron.pfIsolationVariables().sumPhotonEt -
//                                0.5 * patElectron.pfIsolationVariables().sumPUPt, 0.0)) / patElectron.pt();
//
//         syncTree.id_e_mva_nt_loose_1() = Run2::DefaultFillValueForSyncTree();
//
//         syncTree.gen_match_1() = true;
//         // Leg 2, tau
//         const pat::Tau& patTau = selection.GetLeg(2)->GetNtupleObject<pat::Tau>();
//
//         syncTree.pt_2()     = selection.GetLeg(2)->GetMomentum().Pt();
//         syncTree.phi_2()    = selection.GetLeg(2)->GetMomentum().Phi();
//         syncTree.eta_2()    = selection.GetLeg(2)->GetMomentum().Eta();
//         syncTree.m_2()      = selection.GetLeg(2)->GetMomentum().M();
//         syncTree.q_2()      = selection.GetLeg(2)->GetCharge();
//         syncTree.pfmt_2()     = Calculate_MT(selection.GetLeg(2)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
//         syncTree.d0_2()     = Calculate_dxy(selection.GetLeg(2)->GetVertexPosition(), primaryVertex->GetPosition(),
//                                              selection.GetLeg(2)->GetMomentum());
//         syncTree.dZ_2()     = selection.GetLeg(2)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
//         syncTree.iso_2()    = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
//         syncTree.id_e_mva_nt_loose_1() = Run2::DefaultFillValueForSyncTree();
//         syncTree.gen_match_2() = true;
//
//         syncTree.againstElectronLooseMVA6_2()   = patTau.tauID("againstElectronLooseMVA6");
//         syncTree.againstElectronMediumMVA6_2()  = patTau.tauID("againstElectronMediumMVA6");
//         syncTree.againstElectronTightMVA6_2()   = patTau.tauID("againstElectronTightMVA6");
//         syncTree.againstElectronVLooseMVA6_2()  = patTau.tauID("againstElectronVLooseMVA6");
//         syncTree.againstElectronVTightMVA6_2()  = patTau.tauID("againstElectronVTightMVA6");
//
//         syncTree.againstMuonLoose3_2()          = patTau.tauID("againstMuonLoose3");
//         syncTree.againstMuonTight3_2()          = patTau.tauID("againstMuonTight3");
//
//         syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = patTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
//         syncTree.byIsolationMVA3newDMwLTraw_2()               = patTau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
//         syncTree.byIsolationMVA3oldDMwLTraw_2()               = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
//         syncTree.byIsolationMVA3newDMwoLTraw_2()              = Run2::DefaultFillValueForSyncTree();
//         syncTree.byIsolationMVA3oldDMwoLTraw_2()              = Run2::DefaultFillValueForSyncTree();
//
//         syncTree.byVLooseIsolationMVArun2v1DBoldDMwLT_2()     = patTau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
//         syncTree.byLooseIsolationMVArun2v1DBoldDMwLT_2()      = patTau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
//         syncTree.byMediumIsolationMVArun2v1DBoldDMwLT_2()     = patTau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
//         syncTree.byTightIsolationMVArun2v1DBoldDMwLT_2()      = patTau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
//         syncTree.byVTightIsolationMVArun2v1DBoldDMwLT_2()     = patTau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
//
//         syncTree.decayModeFindingOldDMs_2() = patTau.tauID("decayModeFinding");
//
//         syncTree.dilepton_veto()  = selection.Zveto;
//         syncTree.extraelec_veto() = selection.electronVeto;
//         syncTree.extramuon_veto() = selection.muonVeto;
//
//         // Jets
//         syncTree.njetspt20() = selection.jets.size();
//         //syncTree.njets()     = selection.numJet;
//         syncTree.nbtag()     = selection.bjets.size();
//         //int jetCount = 0;
//
//         Int_t numJet = 0;
//         for( const CandidateV2Ptr& jet : selection.jets ){
//                 const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
//                 if (jet->GetMomentum().Pt() > 30 ) numJet++;
//                 syncTree.pt_jets()   .push_back(jet->GetMomentum().Pt());
//                 syncTree.eta_jets()  .push_back(jet->GetMomentum().Eta());
//                 syncTree.phi_jets()  .push_back(jet->GetMomentum().Phi());
//                 syncTree.energy_jets() .push_back(jet->GetMomentum().E());
//                 syncTree.rawf_jets() .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
//                 syncTree.mva_jets()  .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
//                 syncTree.csv_jets()  .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
//                 syncTree.partonFlavour_jets() .push_back(pat_jet1.partonFlavour());
//          }
//         syncTree.njets() = numJet;
//
//
//         for( const CandidateV2Ptr& jet : selection.bjets ){
//                 const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
//                 syncTree.pt_bjets()   .push_back(jet->GetMomentum().Pt());
//                 syncTree.eta_bjets()  .push_back(jet->GetMomentum().Eta());
//                 syncTree.phi_bjets()  .push_back(jet->GetMomentum().Phi());
//                 syncTree.energy_bjets() .push_back(jet->GetMomentum().E());
//                 syncTree.rawf_bjets() .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
//                 syncTree.mva_bjets()  .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
//                 syncTree.csv_bjets()  .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
//                 syncTree.partonFlavour_bjets() .push_back(pat_jet1.partonFlavour());
//          }
//
//         syncTree.Fill();
    }


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SyncTreeProducer_tautau);
