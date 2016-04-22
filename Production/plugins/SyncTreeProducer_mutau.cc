/*!
 * \file SyncTreeProducer_mutau.cc
 * \author author: Claudio Caputo
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "h-tautau/Production/interface/SyncTreeProducer_mutau.h"

//
// constructors and destructor
//
SyncTreeProducer_mutau::SyncTreeProducer_mutau(const edm::ParameterSet& iConfig):
  BaseEDAnalyzer(iConfig),
  syncTree(&edm::Service<TFileService>()->file(),false),
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
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::MuTau;

  const auto Key = analysis::stringToDataSourceTypeMap.at(BaseEDAnalyzer::GetSampleType());
  const auto& hltPaths = MuTau::trigger::hltPathMaps.at(Key);

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

      //Di-Lepton Veto
	  const auto z_muons = CollectZmuons();
	  const auto z_muons_candidates = BaseEDAnalyzer::FindCompatibleObjects(z_muons,z_muons,DeltaR_betweenSignalObjects,analysis::CandidateV2::Type::Z,"Z_mu_mu",0);
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

	  analysis::MissingETPtr pfMET(new analysis::MissingET((*pfMETs)[0],*(BaseEDAnalyzer::GetMETCovMatrix())));
      selection.pfMET = pfMET;

      selection.jets   = CollectJets();
      selection.bjets  = CollectBJets();

//      double svfit_mass = SVFit(higgs,pfMET);
//          std::cout<< "\n\t CHOOSEN SVFit -->  " << svfit_mass <<std::endl;
      selection.svfitResult = BaseEDAnalyzer::SVFit<pat::Muon>(higgs,pfMET);
          std::cout<< "\n\t CHOOSEN SVFit -->  " << selection.svfitResult.mass <<std::endl;


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

void SyncTreeProducer_mutau::SelectJets(const analysis::CandidateV2Ptr& jet, analysis::SelectionManager& selectionManager, cuts::Cutter& cut)
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

void SyncTreeProducer_mutau::FillSyncTree(const edm::Event& iEvent)
    {
    	using namespace analysis;
        static const float default_value = Run2::DefaultFloatFillValueForSyncTree();


		    const VertexV2Ptr primaryVertex(new VertexV2((*(BaseEDAnalyzer::GetVertexCollection())).front()));

        //BaseEDAnalyzer::FillSyncTree(iEvent, selection, syncTree);

        //const auto primaryVertex = (*(BaseEDAnalyzer::GetVertexCollection())).ptrAt(0);

        // Event
        std::cout<<"~~~~~~~~~~~~~EVENT Info~~~~~~~~~"<<std::endl;
        std::cout<<"Run = "<<iEvent.id().run()<<"  Lumi = "<<iEvent.id().luminosityBlock()
                <<" Event = "<<iEvent.id().event()<<std::endl;
        syncTree().run  = iEvent.id().run();
        syncTree().lumi = iEvent.id().luminosityBlock();
        syncTree().evt  = iEvent.id().event();
        syncTree().channelID = static_cast<int>(SyncTreeProducer_mutau::ChannelId());


        syncTree().HT   = selection.HT;
        if(selection.HT<100) syncTree().HTBin = static_cast<int>(Run2::HTbinning::lt100);
        if(100<=selection.HT && selection.HT<200) syncTree().HTBin = static_cast<int>(Run2::HTbinning::f100to200);
        if(200<=selection.HT && selection.HT<400) syncTree().HTBin = static_cast<int>(Run2::HTbinning::f200to400);
        if(400<=selection.HT && selection.HT<600) syncTree().HTBin = static_cast<int>(Run2::HTbinning::f400to600);
        if(selection.HT>=600) syncTree().HTBin = static_cast<int>(Run2::HTbinning::gt600);
        if(selection.HT<0)    syncTree().HTBin = -1;

         if (BaseEDAnalyzer::isMC()) syncTree().weightevt = selection.weightevt;
         else syncTree().weightevt = default_value;

        syncTree().npv = selection.npv;
        syncTree().npu = selection.numtruepileupinteractions;


        // HTT candidate
        syncTree().m_vis = selection.higgs->GetMomentum().M();
        syncTree().pt_tt  = (selection.GetLeg(1)->GetMomentum() + selection.GetLeg(2)->GetMomentum()).Pt();
        syncTree().m_sv = selection.svfitResult.has_valid_mass
                ? selection.svfitResult.mass : Run2::DefaultFillValueForSyncTree();
        syncTree().pt_sv = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Pt() : Run2::DefaultFillValueForSyncTree();
        syncTree().eta_sv = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Eta() : Run2::DefaultFillValueForSyncTree();
        syncTree().phi_sv = selection.svfitResult.has_valid_momentum
                ? selection.svfitResult.momentum.Phi() : Run2::DefaultFillValueForSyncTree();


        syncTree().met        = selection.pfMET->Pt();
        syncTree().metphi     = selection.pfMET->Phi();
        syncTree().isPFMET    = selection.pfMET->isPFMET();
        syncTree().metcov00   = selection.pfMET->GetCovVector().at(0);
        syncTree().metcov01   = selection.pfMET->GetCovVector().at(1);
        syncTree().metcov10   = selection.pfMET->GetCovVector().at(2);
        syncTree().metcov11   = selection.pfMET->GetCovVector().at(3);

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
        syncTree().gen_match_1 = true;
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
        syncTree().gen_match_2 = true;

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

        syncTree().dilepton_veto  = selection.Zveto;
        syncTree().extraelec_veto = selection.electronVeto;
        syncTree().extramuon_veto = selection.muonVeto;

        // Jets
        syncTree().njetspt20 = selection.jets.size();
        //syncTree.njets()     = selection.numJet;
        syncTree().nbtag     = selection.bjets.size();
        //int jetCount = 0;

        Int_t numJet = 0;
        for( const CandidateV2Ptr& jet : selection.jets ){
                const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
                if (jet->GetMomentum().Pt() > 30 ) numJet++;
                syncTree().pt_jets     .push_back(jet->GetMomentum().Pt());
                syncTree().eta_jets    .push_back(jet->GetMomentum().Eta());
                syncTree().phi_jets    .push_back(jet->GetMomentum().Phi());
                syncTree().energy_jets .push_back(jet->GetMomentum().E());
                syncTree().rawf_jets   .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
                syncTree().mva_jets    .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
                syncTree().csv_jets    .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                syncTree().partonFlavour_jets .push_back(pat_jet1.partonFlavour());
         }
        syncTree().njets = numJet;


        for( const CandidateV2Ptr& jet : selection.bjets ){
                const pat::Jet& pat_jet1 = jet->GetNtupleObject<pat::Jet>();
                syncTree().pt_bjets     .push_back(jet->GetMomentum().Pt());
                syncTree().eta_bjets    .push_back(jet->GetMomentum().Eta());
                syncTree().phi_bjets    .push_back(jet->GetMomentum().Phi());
                syncTree().energy_bjets .push_back(jet->GetMomentum().E());
                syncTree().rawf_bjets   .push_back((pat_jet1.correctedJet("Uncorrected").pt() ) / jet->GetMomentum().Pt());
                syncTree().mva_bjets    .push_back(pat_jet1.userFloat("pileupJetId:fullDiscriminant"));
                syncTree().csv_bjets    .push_back(pat_jet1.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                syncTree().partonFlavour_bjets .push_back(pat_jet1.partonFlavour());
         }

         syncTree.Fill();
    }

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SyncTreeProducer_mutau);
