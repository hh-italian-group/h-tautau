/*! Implementation of a SyncTree producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Production/interface/SyncTreeProducer_mutau.h"

//
// constructors and destructor
//
SyncTreeProducer_mutau::SyncTreeProducer_mutau(const edm::ParameterSet& iConfig):
  BaseEDAnalyzer(iConfig),
  // syncTree(&edm::Service<TFileService>()->file(),false),
  sync_tree(new ntuple::SyncTree("mutau",&edm::Service<TFileService>()->file(),false)),
  syncTree(*sync_tree),
  anaData("mutau_cuts.root"){}


SyncTreeProducer_mutau::~SyncTreeProducer_mutau()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

analysis::Channel SyncTreeProducer_mutau::ChannelId() const {return analysis::Channel::MuTau;}
//

// member functions
//

// ------------ method called for each event  ------------
void
SyncTreeProducer_mutau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::Central);
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::TauUp);
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::TauDown);
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::JetUp);
  ProcessEvent(iEvent,iSetup,analysis::EventEnergyScale::JetDown);

}

void
SyncTreeProducer_mutau::ProcessEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, const analysis::EventEnergyScale eventEnergyScale)
{
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::MuTau;

  const auto Key = analysis::stringToDataSourceTypeMap.at(BaseEDAnalyzer::GetSampleType());
  const auto& hltPaths = MuTau::trigger::hltPathMaps.at(Key);
 
  /*  Jet Scale Up and Down   -  to be moved in BaseEDAnalyzer*/
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get("AK5PF",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  jecUnc = new JetCorrectionUncertainty(JetCorPar);
  /*--------------------------*/ 
   
  cuts::Cutter cut(&GetAnaData().Selection("events"));

  //Get collection
  Initialize(iEvent);

  try{

      selection.eventEnergyScale = eventEnergyScale;
      
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
	  const auto z_muons = CollectZmuons();
	  const auto z_muons_candidates = BaseEDAnalyzer::FindCompatibleObjects(z_muons,z_muons,DeltaR_DileptonVeto,analysis::CandidateV2::Type::Z,"Z_mu_mu",0);
      selection.Zveto = z_muons_candidates.size() ? true : false;

      //Signal-like leptons selection
      const auto selectedMuons = CollectSignalMuons();

      cut(selectedMuons.size(),"muons");
      //std::cout<<" Vertici    --->   "<<PV->ndof()<<std::endl;

      const auto selectedTaus = CollectSignalTaus();

      cut(selectedTaus.size(),"taus");

      auto higgses = BaseEDAnalyzer::FindCompatibleObjects(selectedMuons,selectedTaus,DeltaR_betweenSignalObjects,analysis::CandidateV2::Type::Higgs,
                                         "H_mu_tau");

      cut(higgses.size(),"mu_tau_pair");

      bool isCrossTrigger = false;
      auto triggeredHiggses = ApplyTriggerMatch(iEvent, higgses,hltPaths,false,isCrossTrigger);

      cut(triggeredHiggses.size(),"triggerMatch");

      selection.higgs = SelectHiggs(triggeredHiggses);
   	  auto higgs = selection.higgs;

      //Third-Lepton Veto
      const auto muonVetoCollection     = CollectVetoMuons();
      const auto electronVetoCollection = CollectVetoElectrons();

      selection.electronVeto = electronVetoCollection.size() ? true : false;

      selection.muonVeto = false;
      for (auto &muon: muonVetoCollection){
           if (!(selection.GetLeg(1)->GetMomentum().DeltaR(muon->GetMomentum()) < 0.05)) selection.muonVeto = true;
      }

	  
	  /*  MET Shifted by JEC */
	  analysis::MissingETPtr pfMET(new analysis::MissingET((*pfMETs)[0],*(BaseEDAnalyzer::GetMETCovMatrix())));
     if ( selection.eventEnergyScale == analysis::EventEnergyScale::JetUp || selection.eventEnergyScale == analysis::EventEnergyScale::JetDown ){
         std::map<analysis::EventEnergyScale, pat::MET::METUncertainty > metUncertantyMap = {{analysis::EventEnergyScale::JetUp,pat::MET::METUncertainty::JetEnUp},
                                                                                             {analysis::EventEnergyScale::JetDown,pat::MET::METUncertainty::JetEnDown}};
         const reco::Candidate::LorentzVector shifted_p4 = ((*pfMETs)[0]).shiftedP4(metUncertantyMap[selection.eventEnergyScale]);
         std::cout<<" MET   P4 Value:  "<< shifted_p4 <<std::endl;
         pfMET->ShiftMET(shifted_p4);
     }
      selection.pfMET = pfMET;

      selection.jets        = CollectJets();
      //selection.jetsTight   = CollectJetsTight();
      selection.bjets       = CollectBJets();
      
      // Sorting for KinFit
      // -----------
      // auto jetCsvOrdering = [&](const analysis::CandidateV2Ptr first,const analysis::CandidateV2Ptr second) -> bool { const pat::Jet& pat_jet1 = first->GetNtupleObject<pat::Jet>();
//                                                                                                                       const pat::Jet& pat_jet2 = second->GetNtupleObject<pat::Jet>();
//                                                                                                                       return pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > pat_jet2.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");};
//       std::sort(selection.jets.begin(),selection.jets.end(),jetCsvOrdering);
//       std::sort(selection.jetsTight.begin(),selection.jetsTight.end(),jetCsvOrdering);
//       std::sort(selection.bjets.begin(),selection.bjets.end(),jetCsvOrdering);

      selection.svfitResult = BaseEDAnalyzer::SVFit<pat::Muon>(higgs,pfMET);
      std::cout<< "\n\t CHOOSEN SVFit -->  " << selection.svfitResult.mass <<std::endl;

     //if (selection.jetsTight.size()>=2) GetKinFitResults();

	  FillSyncTree(iEvent);

  }catch(cuts::cut_failed&){}

  GetAnaData().Selection("events").fill_selection();

}

// ------------ objects selection  ------------
void SyncTreeProducer_mutau::SelectSignalMuon(const analysis::CandidateV2Ptr& muon, analysis::SelectionManager& selectionManager, cuts::Cutter& cut){
	using namespace cuts::Htautau_2015::MuTau;
    using namespace cuts::Htautau_2015::MuTau::muonID;
    const pat::Muon& object = muon->GetNtupleObject<pat::Muon>();
    const auto PV = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);

    // bool muonIP = fabs(object.muonBestTrack()->dxy(PV->position())) < muonID::dB &&
//                   fabs(object.muonBestTrack()->dz(PV->position())) < muonID::dz;

    cut(true, ">0 mu cand");
    cut(X(pt()) > pt, "pt");
    cut(std::abs( X(eta()) ) < eta, "eta");
    const double muonDB = fabs(object.muonBestTrack()->dxy(PV->position()));
    cut(Y(muonDB) < muonID::dB, "dxy");
    const double muonDZ = fabs(object.muonBestTrack()->dz(PV->position()));
    cut(Y(muonDZ) < muonID::dz, "dz");
    cut(X(isMediumMuon()) == muonID::isMediumMuon, "muonID");
}

void SyncTreeProducer_mutau::SelectSignalTau(const analysis::CandidateV2Ptr& tau, analysis::SelectionManager& selectionManager, cuts::Cutter& cut){
	using namespace cuts::Htautau_2015::MuTau;
    using namespace cuts::Htautau_2015::MuTau::tauID;
    const pat::Tau& object = tau->GetNtupleObject<pat::Tau>();
    
    //std::cout<<" @@@@@@@@  Tau pt before ES : " << tau->GetMomentum() << std::endl;
    
    if(selection.eventEnergyScale == analysis::EventEnergyScale::TauUp || selection.eventEnergyScale == analysis::EventEnergyScale::TauDown) {
            const double sign = selection.eventEnergyScale == analysis::EventEnergyScale::TauUp ? +1 : -1;
            const double sf = 1.0 + sign * cuts::Htautau_2015::tauCorrections::energyUncertainty;
            tau->ScaleMomentum(sf);
    }
    //    std::cout<<" @@@@@@@@  Tau pt after ES : " << tau->GetMomentum() << std::endl;

    TLorentzVector momentum= tau->GetMomentum();
    // if(selection.eventEnergyScale == analysis::EventEnergyScale::TauUp)   momentum = tau->GetScaledUpMomentum();
//     if(selection.eventEnergyScale == analysis::EventEnergyScale::TauDown) momentum = tau->GetScaledDownMomentum();
    
    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(object.leadChargedHadrCand().get());

    double const taudz = packedLeadTauCand->dz();
    int const decayModeFinding    = object.tauID("decayModeFinding");

    cut(true, ">0 tau cand");
//     cut(X(pt()) > pt, "pt");
//     cut(std::fabs( X(eta()) ) < eta, "eta");
    cut(Y(momentum.Pt()) > pt, "pt");
    cut(std::fabs( Y(momentum.Eta()) ) < eta, "eta");
    cut(Y(decayModeFinding) > tauID::decayModeFinding, "oldDecayMode");
    cut(std::fabs( Y(taudz) ) < dz ,"dz");
    cut(std::abs(X(charge()))==1,"charge");
}

void SyncTreeProducer_mutau::SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
    {
        using namespace cuts::Htautau_2015;
        using namespace cuts::Htautau_2015::jetID;
        const pat::Jet& object = jet->GetNtupleObject<pat::Jet>();
        
        if(selection.eventEnergyScale == analysis::EventEnergyScale::JetUp || selection.eventEnergyScale == analysis::EventEnergyScale::JetUp) {
            jecUnc->setJetEta(X(eta()));
            jecUnc->setJetPt(X(pt())); // here you must use the CORRECTED jet pt
            double unc = jecUnc->getUncertainty(true);
            const double sign = selection.eventEnergyScale == analysis::EventEnergyScale::JetUp ? +1 : -1;
            const double sf = (1+(sign*unc)) ; // shift = +1(up), or -1(down)
            jet->ScaleMomentum(sf);
        }
         
         TLorentzVector momentum= jet->GetMomentum();
        //std::cout<<" \t\t\t\t Jet Pt  :   "<< X(pt()) <<"  Jet Pt Shifted  : "<<ptCor_shifted<<std::endl;
 
        cut(true, ">0 jet cand");
        cut(Y(momentum.Pt()) > pt_loose, "pt_loose");
        cut(std::abs( Y(momentum.Eta()) ) < eta, "eta");
        const bool jetPFID = passPFLooseId(object);
        cut(Y(jetPFID),"jet_id");
        const double deltaR_leg1 = jet->GetMomentum().DeltaR(selection.GetLeg(1)->GetMomentum());
        cut(Y(deltaR_leg1) > deltaR_signalObjects,"deltaR_lep1");
        const double deltaR_leg2 = jet->GetMomentum().DeltaR(selection.GetLeg(2)->GetMomentum());
        cut(Y(deltaR_leg2) > deltaR_signalObjects,"deltaR_lep2");
    }

void SyncTreeProducer_mutau::SelectJetsTight(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
{
   using namespace cuts::Htautau_2015;
    using namespace cuts::Htautau_2015::jetID;
    //const pat::Jet& object = jet->GetNtupleObject<pat::Jet>();
   
    SelectJets(jet,selectionManager, cut);
    TLorentzVector momentum= jet->GetMomentum();
    cut(std::abs( Y(momentum.Pt()) ) > 30, "eta");
    cut(std::abs( Y(momentum.Eta()) ) < btag::eta, "eta");
}
void SyncTreeProducer_mutau::SelectBJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
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
SyncTreeProducer_mutau::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SyncTreeProducer_mutau::endJob()
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
// void
// SyncTreeProducer_mutau::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }


analysis::CandidateV2Ptr SyncTreeProducer_mutau::SelectHiggs(analysis::CandidateV2PtrVector& higgses)
{
    if(!higgses.size())
        throw std::runtime_error("no available higgs candidate to select");
    if(higgses.size()==1) return higgses.front();

    const auto higgsSelector = [&] (const analysis::CandidateV2Ptr& first, const analysis::CandidateV2Ptr& second) -> bool
    {
        const pat::Muon& first_muon = first->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
        const pat::Tau&  first_tau  = first->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();
        const pat::Muon& second_muon = second->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
        const pat::Tau&  second_tau  = second->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();

        double iso_mu1 = (first_muon.pfIsolationR03().sumChargedHadronPt + std::max(
                          first_muon.pfIsolationR03().sumNeutralHadronEt +
                          first_muon.pfIsolationR03().sumPhotonEt -
                          0.5 * first_muon.pfIsolationR03().sumPUPt, 0.0)) / first_muon.pt();
        double iso_mu2 = (second_muon.pfIsolationR03().sumChargedHadronPt + std::max(
                          second_muon.pfIsolationR03().sumNeutralHadronEt +
                          second_muon.pfIsolationR03().sumPhotonEt -
                          0.5 * second_muon.pfIsolationR03().sumPUPt, 0.0)) / second_muon.pt();

        double iso_tau1 = first_tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
        double iso_tau2 = second_tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");

        std::cout << "LAMBDA ------------------------  \n IsoMu1 = "<<iso_mu1<<"  IsoMu2 = "<<iso_mu2
                  << "  IsoTau1 = "<<iso_tau1<<"  IsoTau2 = "<<iso_tau2<<std::endl;
        std::cout << " PtMu1 = "<<first_muon.pt()<<"  PtMu2 = "<<second_muon.pt()
                  << "  PtTau1 = "<<first_tau.pt()<<"  PtTau2 = "<<second_tau.pt()<<std::endl;

        if (!iso_mu1 || !iso_mu2 ){
            std::cout << "\t 1st Muon \t ChargedHadronPt = "<<first_muon.pfIsolationR03().sumChargedHadronPt
                      << "\t NeutralHadronEt = "<<first_muon.pfIsolationR03().sumNeutralHadronEt
                      << "\t PhotonEt = "<<first_muon.pfIsolationR03().sumPhotonEt
                      << "\t PUPt = "<<first_muon.pfIsolationR03().sumPUPt<<std::endl;
            std::cout << "\t 2nd Muon \t ChargedHadronPt = "<<second_muon.pfIsolationR03().sumChargedHadronPt
                      << "\t NeutralHadronEt = "<<second_muon.pfIsolationR03().sumNeutralHadronEt
                      << "\t PhotonEt = "<<second_muon.pfIsolationR03().sumPhotonEt
                      << "\t PUPt = "<<second_muon.pfIsolationR03().sumPUPt<<std::endl;
        }
        bool muon = (iso_mu1 < iso_mu2) ||
                    ((iso_mu1 == iso_mu2) ? first_muon.pt() > second_muon.pt() : false);
        bool tau = (iso_tau1 > iso_tau2) ||
                    ((iso_tau1 == iso_tau2) ? first_tau.pt() > second_tau.pt() : false);

//        if ( &first_muon==&second_muon ){
//            if ( &first_tau==&second_tau ) return true;
//            return tau;
//        }

//        if (!muon) return tau;
//        return muon;
        if ( &first_muon!=&second_muon ) return muon;
        if ( &first_tau!=&second_tau ) return tau;

        throw analysis::exception("not found a good criteria for best tau pair");
    };

    std::sort(higgses.begin(), higgses.end(), higgsSelector) ;
    return higgses.front();
}

void SyncTreeProducer_mutau::GetKinFitResults()
    {
        //if(!selection.kinFit_result) {
            std::cout<<"   KinFIT --------------,,,,,,"<<std::endl;
            selection.kinFit_result = std::shared_ptr<analysis::KinFitResult>(new analysis::KinFitResult());
            if(selection.jetsTight.size()>=2) {
               std::cout<<"   ------- KinFIT --------------,,,,,,"<<std::endl;
            
                // if(verbosity > 1) {
//                     std::cout << "b1: " << bjet_momentums.at(selected_bjets.first)
//                               << "\nb2:" << bjet_momentums.at(selected_bjets.second)
//                               << "\ntau1:" << lepton_momentums.at(0)
//                               << "\ntau2:" << lepton_momentums.at(1)
//                               << "\nMET: (" << MET.X() << ", " << MET.Y() << ")"
//                               << "\nMET cov:" << MET_covariance << std::endl;
//                 }
               TMatrixD covMET(2, 2);
               covMET[0][0] = selection.pfMET->GetCovVector().at(0);
               covMET[1][0] = selection.pfMET->GetCovVector().at(1);
               covMET[0][1] = selection.pfMET->GetCovVector().at(2);
               covMET[1][1] = selection.pfMET->GetCovVector().at(3);

                HHKinFit2::HHKinFitMasterHeavyHiggs kin_fit(selection.jetsTight.at(0)->GetMomentum(),
                                                            selection.jetsTight.at(1)->GetMomentum(),
                                                            selection.GetLeg(1)->GetMomentum(), 
                                                            selection.GetLeg(2)->GetMomentum(),
                                                            TVector2(selection.pfMET->Px(), selection.pfMET->Py()), 
                                                            covMET);
                kin_fit.verbosity = 1;
                kin_fit.fit();
                selection.kinFit_result->convergence = kin_fit.getConvergence();
                if(selection.kinFit_result->HasMass()) {
                    selection.kinFit_result->mass = kin_fit.getMH();
                    selection.kinFit_result->chi2 = kin_fit.getChi2();
                    selection.kinFit_result->probability = kin_fit.getFitProb();
                    std::cout<<"   ------- KinFIT --------------,,,,,,"<< selection.kinFit_result->mass <<std::endl;

                }
            }
       // }
    }

void SyncTreeProducer_mutau::FillSyncTree(const edm::Event& iEvent)
    {
         using namespace analysis;
        //static const float default_value = Run2::DefaultFloatFillValueForSyncTree();


		    const VertexV2Ptr primaryVertex(new VertexV2((*(BaseEDAnalyzer::GetVertexCollection())).front()));
         
         BaseEDAnalyzer::FillSyncTree(iEvent,selection,syncTree);
         

        // Leg 1, lepton
        const pat::Muon& patMuon = selection.GetLeg(1)->GetNtupleObject<pat::Muon>();

        syncTree().pt_1     = selection.GetLeg(1)->GetMomentum().Pt();
        syncTree().phi_1    = selection.GetLeg(1)->GetMomentum().Phi();
        syncTree().eta_1    = selection.GetLeg(1)->GetMomentum().Eta();
        syncTree().m_1      = selection.GetLeg(1)->GetMomentum().M();
        syncTree().q_1      = selection.GetLeg(1)->GetCharge();
        syncTree().pfmt_1   = Calculate_MT(selection.GetLeg(1)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
        syncTree().d0_1     = Calculate_dxy(selection.GetLeg(1)->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg(1)->GetMomentum());
        syncTree().dZ_1     = selection.GetLeg(1)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
        syncTree().iso_1    = (patMuon.pfIsolationR03().sumChargedHadronPt + std::max(
                                   patMuon.pfIsolationR03().sumNeutralHadronEt +
                                   patMuon.pfIsolationR03().sumPhotonEt -
                                   0.5 * patMuon.pfIsolationR03().sumPUPt, 0.0)) / patMuon.pt();

        syncTree().id_e_mva_nt_loose_1 = Run2::DefaultFillValueForSyncTree();
        if (BaseEDAnalyzer::isMC()) 
               syncTree().gen_match_1 = analysis::genMatch(patMuon.p4(),*(BaseEDAnalyzer::GetGenParticles()));
        // Leg 2, tau
        const pat::Tau& patTau = selection.GetLeg(2)->GetNtupleObject<pat::Tau>();

        syncTree().pt_2     = selection.GetLeg(2)->GetMomentum().Pt();
        syncTree().phi_2    = selection.GetLeg(2)->GetMomentum().Phi();
        syncTree().eta_2    = selection.GetLeg(2)->GetMomentum().Eta();
        syncTree().m_2      = selection.GetLeg(2)->GetMomentum().M();
        syncTree().q_2      = selection.GetLeg(2)->GetCharge();
        syncTree().pfmt_2   = Calculate_MT(selection.GetLeg(2)->GetMomentum(), selection.pfMET->Pt(), selection.pfMET->Phi());
        syncTree().d0_2     = Calculate_dxy(selection.GetLeg(2)->GetVertexPosition(), primaryVertex->GetPosition(),
                                             selection.GetLeg(2)->GetMomentum());
        syncTree().dZ_2       = selection.GetLeg(2)->GetVertexPosition().Z() - primaryVertex->GetPosition().Z();
		    syncTree().iso_2      = patTau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");

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
DEFINE_FWK_MODULE(SyncTreeProducer_mutau);
