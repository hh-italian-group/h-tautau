/*! Implementation of a SyncTree producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Production/interface/SyncTreeProducer_etau.h"

//
// constructors and destructor
//
SyncTreeProducer_etau::SyncTreeProducer_etau(const edm::ParameterSet& iConfig):
  BaseEDAnalyzer(iConfig),
  //sync_tree(BaseEDAnalyzer::GetSyncTree(&edm::Service<TFileService>()->file())),
  sync_tree(new ntuple::SyncTree("etau",&edm::Service<TFileService>()->file(),false)),
  syncTree(*sync_tree),
  anaData("etau_cuts.root") {}


SyncTreeProducer_etau::~SyncTreeProducer_etau()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

analysis::Channel SyncTreeProducer_etau::ChannelId() const {return analysis::Channel::ETau;}

//
// member functions
//

// ------------ method called for each event  ------------
void
SyncTreeProducer_etau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::Central);
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::TauUp);
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::TauDown);

}

void
SyncTreeProducer_etau::ProcessEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, const analysis::EventEnergyScale tauEnergyScale)
{
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::ETau;

  const auto Key = analysis::stringToDataSourceTypeMap.at(BaseEDAnalyzer::GetSampleType());
  const auto& hltPaths = ETau::trigger::hltPathMaps.at(Key);

  cuts::Cutter cut(&GetAnaData().Selection("events"));

  //Get collection
  Initialize(iEvent);

  try{

      selection.tauEnergyScale = tauEnergyScale;
            
      cut(true,"events");

      cut(HaveTriggerFired(iEvent, hltPaths),"trigger");

      if(BaseEDAnalyzer::isMC()){
         auto lheInfoPair    = BaseEDAnalyzer::computeHtValue();
         selection.HT        = lheInfoPair.first;
         selection.NOutPartons = lheInfoPair.second;
         selection.weightevt = (BaseEDAnalyzer::GetGenEventInfo())->weight();
      }

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

      //Di-Lepton Veto
	  const auto z_electrons = CollectZelectrons();
	  const auto z_electrons_candidates = BaseEDAnalyzer::FindCompatibleObjects(z_electrons,z_electrons,DeltaR_DileptonVeto,analysis::CandidateV2::Type::Z,"Z_mu_mu",0);
     selection.Zveto = z_electrons_candidates.size() ? true : false;

      const auto selectedElectrons = CollectSignalElectrons();

      cut(selectedElectrons.size(),"electrons");

      const auto selectedTaus = CollectSignalTaus();

      cut(selectedTaus.size(),"taus");

      auto higgses = BaseEDAnalyzer::FindCompatibleObjects(selectedElectrons,selectedTaus,DeltaR_betweenSignalObjects,analysis::CandidateV2::Type::Higgs,
                                         "H_mu_tau");

      cut(higgses.size(),"ele_tau_pair");

      bool isCrossTrigger = false;
      auto triggeredHiggses = ApplyTriggerMatch(iEvent, higgses,hltPaths,false,isCrossTrigger);

      cut(triggeredHiggses.size(),"triggerMatch");

      selection.higgs = SelectHiggs(triggeredHiggses);
      auto higgs = selection.higgs;

      //Third-Lepton Veto
      const auto muonVetoCollection     = CollectVetoMuons();
      const auto electronVetoCollection = CollectVetoElectrons();

      selection.muonVeto = muonVetoCollection.size() ? true : false;

      selection.electronVeto = false;
      for (auto &electron: electronVetoCollection){
           if (!(selection.GetLeg(1)->GetMomentum().DeltaR(electron->GetMomentum()) < 0.05)) selection.electronVeto = true;
      }

      analysis::MissingETPtr pfMET(new analysis::MissingET((*pfMETs)[0],*(BaseEDAnalyzer::GetMETCovMatrix())));
      selection.pfMET = pfMET;

      selection.jets   = CollectJets();
      selection.bjets  = CollectBJets();

      selection.svfitResult = BaseEDAnalyzer::SVFit<pat::Electron>(higgs,pfMET);
          std::cout<< "\n\t CHOOSEN SVFit -->  " << selection.svfitResult.mass <<std::endl;

      FillSyncTree(iEvent);

  }catch(cuts::cut_failed&){}

  GetAnaData().Selection("events").fill_selection();

}

// ------------ objects selection  ------------
void SyncTreeProducer_etau::SelectSignalElectron(const analysis::CandidateV2Ptr& electron, analysis::SelectionManager& selectionManager, cuts::Cutter& cut){
    using namespace cuts::Htautau_2015::ETau;
    using namespace cuts::Htautau_2015::ETau::electronID;
    const pat::Electron& object = electron->GetNtupleObject<pat::Electron>();

    const auto PV = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);
    const Ptr<pat::Electron> elPtr(electrons, electron->GetIndex() );


    cut(true, ">0 mu cand");
//    std::cout<<"  edm::Ptr  id  =  "<<elPtr.id()<<"  edm::Ptr  key  =  "<<elPtr.key()<<
//               " MVA TightID  = "<<isTight<<std::endl;

    cut(X(pt()) > pt, "pt");
    cut(std::abs( X(eta()) ) < eta, "eta");
    const double electronD0 = fabs(object.gsfTrack()->dxy(PV->position()));
    cut(Y(electronD0) < electronID::d0, "dxy");
    const double electronDZ = fabs(object.gsfTrack()->dz(PV->position()));
    cut(Y(electronDZ) < electronID::dz, "dz");
    const bool isTight  = (*(BaseEDAnalyzer::GetTElectronIDValueMap()))[elPtr];
    cut(isTight, "electronMVATightID");
    const int eleMissingHits = object.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    cut(eleMissingHits <= electronID::missingHits,"missingHits");
    cut(X(passConversionVeto()) == electronID::passConversionVeto,"conversionVeto");
}

void SyncTreeProducer_etau::SelectSignalTau(const analysis::CandidateV2Ptr& tau, analysis::SelectionManager& selectionManager, cuts::Cutter& cut){
    using namespace cuts::Htautau_2015::ETau;
    using namespace cuts::Htautau_2015::ETau::tauID;
    const pat::Tau& object = tau->GetNtupleObject<pat::Tau>();

    std::cout<<" @@@@@@@@  Tau pt before ES : " << tau->GetMomentum() << std::endl;
    
    if(selection.tauEnergyScale == analysis::EventEnergyScale::TauUp || selection.tauEnergyScale == analysis::EventEnergyScale::TauDown) {
            const double sign = selection.tauEnergyScale == analysis::EventEnergyScale::TauUp ? +1 : -1;
            const double sf = 1.0 + sign * cuts::Htautau_2015::tauCorrections::energyUncertainty;
            tau->ScaleMomentum(sf);
    }
        std::cout<<" @@@@@@@@  Tau pt after ES : " << tau->GetMomentum() << std::endl;

    TLorentzVector momentum= tau->GetMomentum();

    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(object.leadChargedHadrCand().get());

    double const taudz = packedLeadTauCand->dz();
    int const decayModeFinding    = object.tauID("decayModeFinding");

    cut(true, ">0 tau cand");
    cut(Y(momentum.Pt()) > pt, "pt");
    cut(std::fabs( Y(momentum.Eta()) ) < eta, "eta");
    cut(Y(decayModeFinding) > tauID::decayModeFinding, "oldDecayMode");
    cut(std::fabs( Y(taudz) ) < dz ,"dz");
    cut(std::abs(X(charge()))==1,"charge");
}

void SyncTreeProducer_etau::SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
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

void SyncTreeProducer_etau::SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
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
SyncTreeProducer_etau::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SyncTreeProducer_etau::endJob()
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


analysis::CandidateV2Ptr SyncTreeProducer_etau::SelectHiggs(analysis::CandidateV2PtrVector& higgses)
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

void SyncTreeProducer_etau::FillSyncTree(const edm::Event& iEvent)
    {
      using namespace analysis;
        //static const float default_value = Run2::DefaultFloatFillValueForSyncTree();


		    const VertexV2Ptr primaryVertex(new VertexV2((*(BaseEDAnalyzer::GetVertexCollection())).front()));
         
         BaseEDAnalyzer::FillSyncTree(iEvent,selection,syncTree);

        // Leg 1, lepton
        const pat::Electron& patElectron = selection.GetLeg(1)->GetNtupleObject<pat::Electron>();

        syncTree().pt_1     = selection.GetLeg(1)->GetMomentum().Pt();
        syncTree().phi_1    = selection.GetLeg(1)->GetMomentum().Phi();
        syncTree().eta_1    = selection.GetLeg(1)->GetMomentum().Eta();
        syncTree().m_1      = selection.GetLeg(1)->GetMomentum().M();
        syncTree().q_1      = selection.GetLeg(1)->GetCharge();
        syncTree().pfmt_1     = Calculate_MT(selection.GetLeg(1)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
        syncTree().d0_1     = Calculate_dxy(selection.GetLeg(1)->GetVertexPosition(), primaryVertex->GetPosition(),
                                              selection.GetLeg(1)->GetMomentum());
        syncTree().dZ_1     = selection.GetLeg(1)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
        syncTree().iso_1    = (patElectron.pfIsolationVariables().sumChargedHadronPt
                                 + std::max(patElectron.pfIsolationVariables().sumNeutralHadronEt
                                 + patElectron.pfIsolationVariables().sumPhotonEt -
                                 0.5 * patElectron.pfIsolationVariables().sumPUPt, 0.0)) / patElectron.pt();

        syncTree().id_e_mva_nt_loose_1 = Run2::DefaultFillValueForSyncTree();

        if (BaseEDAnalyzer::isMC()) 
               syncTree().gen_match_1 = analysis::genMatch(patElectron.p4(),*(BaseEDAnalyzer::GetGenParticles()));
        // Leg 2, tau
        const pat::Tau& patTau = selection.GetLeg(2)->GetNtupleObject<pat::Tau>();

        syncTree().pt_2     = selection.GetLeg(2)->GetMomentum().Pt();
        syncTree().phi_2    = selection.GetLeg(2)->GetMomentum().Phi();
        syncTree().eta_2    = selection.GetLeg(2)->GetMomentum().Eta();
        syncTree().m_2      = selection.GetLeg(2)->GetMomentum().M();
        syncTree().q_2      = selection.GetLeg(2)->GetCharge();
        syncTree().pfmt_2     = Calculate_MT(selection.GetLeg(2)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
        syncTree().d0_2     = Calculate_dxy(selection.GetLeg(2)->GetVertexPosition(), primaryVertex->GetPosition(),
                                              selection.GetLeg(2)->GetMomentum());
        syncTree().dZ_2     = selection.GetLeg(2)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
        syncTree().iso_2    = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        syncTree().id_e_mva_nt_loose_1 = Run2::DefaultFillValueForSyncTree();
        if (BaseEDAnalyzer::isMC()) 
               syncTree().gen_match_2 = analysis::genMatch(patTau.p4(),*(BaseEDAnalyzer::GetGenParticles()));

        syncTree().againstElectronLooseMVA6_2   = patTau.tauID("againstElectronLooseMVA6");
        syncTree().againstElectronMediumMVA6_2  = patTau.tauID("againstElectronMediumMVA6");
        syncTree().againstElectronTightMVA6_2   = patTau.tauID("againstElectronTightMVA6");
        syncTree().againstElectronVLooseMVA6_2  = patTau.tauID("againstElectronVLooseMVA6");
        syncTree().againstElectronVTightMVA6_2  = patTau.tauID("againstElectronVTightMVA6");

        syncTree().againstMuonLoose3_2          = patTau.tauID("againstMuonLoose3");
        syncTree().againstMuonTight3_2          = patTau.tauID("againstMuonTight3");

        syncTree().byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = patTau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        syncTree().byIsolationMVA3newDMwLTraw_2               = patTau.tauID("byIsolationMVArun2v1DBnewDMwLTraw");
        syncTree().byIsolationMVA3oldDMwLTraw_2               = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        syncTree().byIsolationMVA3newDMwoLTraw_2              = Run2::DefaultFillValueForSyncTree();
        syncTree().byIsolationMVA3oldDMwoLTraw_2              = Run2::DefaultFillValueForSyncTree();

        syncTree().byVLooseIsolationMVArun2v1DBoldDMwLT_2     = patTau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
        syncTree().byLooseIsolationMVArun2v1DBoldDMwLT_2      = patTau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
        syncTree().byMediumIsolationMVArun2v1DBoldDMwLT_2     = patTau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
        syncTree().byTightIsolationMVArun2v1DBoldDMwLT_2      = patTau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
        syncTree().byVTightIsolationMVArun2v1DBoldDMwLT_2     = patTau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

        syncTree().decayModeFindingOldDMs_2 = patTau.tauID("decayModeFinding");

        syncTree.Fill();
    }


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SyncTreeProducer_etau);
