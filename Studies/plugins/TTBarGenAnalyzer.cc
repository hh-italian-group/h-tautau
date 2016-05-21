// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"

#include"h-tautau/Studies/interface/TTBarTree.h"

#include "TMath.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//
namespace analysis{
  struct GenState{
    const reco::Candidate *top;
    const reco::Candidate *W;
    std::vector<const reco::Candidate*> W_daugther;
    const reco::Candidate *bjet, *genJetMatched;
  };

   typedef std::vector<GenState> GenStateVector;

}


class TTBarGenAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TTBarGenAnalyzer(const edm::ParameterSet&);
      ~TTBarGenAnalyzer();
      bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);



   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      const reco::Candidate * findDecayParticle (const reco::Candidate* );
      void FillTree(const edm::Event&);

      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<edm::View<reco::GenJet> > genJetToken_;

      analysis::GenStateVector ttbarStateVector;
      ntuple::TTBarTree tree;
};

TTBarGenAnalyzer::TTBarGenAnalyzer(const edm::ParameterSet& iConfig):
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  genJetToken_(consumes< edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJet"))),
  tree(&edm::Service<TFileService>()->file(),false)
{
}


TTBarGenAnalyzer::~TTBarGenAnalyzer()
{
}

//Check recursively if any ancestor of particle is the given one

bool TTBarGenAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
//particle is already the ancestor
        if(ancestor == particle ) return true;

//otherwise loop on mothers, if any and return true if the ancestor is found
        for(size_t i=0;i< particle->numberOfMothers();i++)
        {
                if(isAncestor(ancestor,particle->mother(i))) return true;
        }
//if we did not return yet, then particle and ancestor are not relatives
        return false;
}

const reco::Candidate * TTBarGenAnalyzer::findDecayParticle (const reco::Candidate* particle){

   if (particle->numberOfDaughters() == 1){
      const reco::Candidate *daughter = particle->daughter(0);

      return findDecayParticle(daughter);
   }

   std::cout << "\t\t New Daughters  :  "<< particle->numberOfDaughters()<<std::endl;
   return particle;
}

void
TTBarGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
        using namespace edm;
        using namespace reco;
        using namespace pat;

        // Pruned particles are the one containing "important" stuff
        Handle<edm::View<reco::GenParticle> > pruned;
        iEvent.getByToken(prunedGenToken_,pruned);

        // Packed particles are all the status 1, so usable to remake jets
        // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
        Handle<edm::View<pat::PackedGenParticle> > packed;
        iEvent.getByToken(packedGenToken_,packed);

        Handle<edm::View<reco::GenJet>> genJet;
        iEvent.getByToken(genJetToken_,genJet);

        //let's try to find all status1 originating directly from a B meson decay
        analysis::GenState ttbarState;

        for(size_t i=0; i<pruned->size();i++){
                // if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600)
                if(abs((*pruned)[i].pdgId()) == 6 ){
                        const Candidate * top = &(*pruned)[i];
                        if ((*pruned)[i].status() == 62){
                          ttbarState.top = top;
                        std::cout << "PdgID: " << top->pdgId() << " pt " << top->pt() << " eta: " << top->eta() << " phi: " << top->phi() << std::endl;
                        std::cout << " Number of Daughters:  " << top->numberOfDaughters() << "  Status: " << top->status() << std::endl;
                        for (size_t i=0; i<top->numberOfDaughters(); ++i){
                           const Candidate * daughter = top->daughter(i);
                           std::cout << "\t\t PdgID: " << daughter->pdgId() << " pt " << daughter->pt() << " eta: " << daughter->eta() << " phi: " << daughter->phi() << std::endl;

                           if ( std::abs(daughter->pdgId()) == 24){
                              auto W = findDecayParticle(daughter);
                              ttbarState.W = W;
                              for (size_t i=0; i<W->numberOfDaughters(); ++i){
                                 const Candidate * Wdaughter = W->daughter(i);
                                 ttbarState.W_daugther.push_back(Wdaughter);
                                 std::cout << "\t\t\t PdgID: " << Wdaughter->pdgId() << " pt " << Wdaughter->pt() << " eta: " << Wdaughter->eta() << " phi: " << Wdaughter->phi() << std::endl;
                              }
                           }
                           if ( std::abs(daughter->pdgId()) == 5 ){
                             ttbarState.bjet = daughter;
                              for (size_t i=0; i<genJet->size();i++){
                                 const GenJet * jet = &(*genJet)[i];
                                 if (ROOT::Math::VectorUtil::DeltaR( daughter->p4(), jet->p4()) <0.02 ) ttbarState.genJetMatched = jet;
                                 std::cout << "\t\t\t PdgID: " << jet->pdgId() << " pt " << jet->pt() << " eta: " << jet->eta() << " phi: " << jet->phi() << std::endl;
                              }
                           }
                        }
                       //    std::cout << "  found daugthers: " << std::endl;
                      //   for(size_t j=0; j<packed->size();j++){
// //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
//                                 const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
//                                 if(motherInPrunedCollection != nullptr && isAncestor( top , motherInPrunedCollection)){
//                                         std::cout << "     PdgID: " << (*packed)[j].pdgId() << " pt " << (*packed)[j].pt() << " eta: " << (*packed)[j].eta() << " phi: " << (*packed)[j].phi() << " Status:  " << (*packed)[j].status() <<std::endl;
//                                 }
//
//                         }
                        ttbarStateVector.push_back(ttbarState);

                       }
                }

        }

      FillTree(iEvent);
	ttbarStateVector.clear();
}


// ------------ method called once each job just before starting event loop  ------------
        void
TTBarGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
        void
TTBarGenAnalyzer::endJob()
{
   tree.Write();
}

void TTBarGenAnalyzer::FillTree(const edm::Event& iEvent){

      tree().run  = iEvent.id().run();
      tree().lumi = iEvent.id().luminosityBlock();
      tree().evt  = iEvent.id().event();

      for (auto &ttbarState : ttbarStateVector) {
        tree().pdgid_tops.push_back(ttbarState.top->pdgId());
        tree().pt_tops.push_back(ttbarState.top->pt());
        tree().eta_tops.push_back(ttbarState.top->eta());
        tree().phi_tops.push_back(ttbarState.top->phi());
        tree().mass_tops.push_back(ttbarState.top->mass());
        tree().energy_tops.push_back(ttbarState.top->energy());

        tree().pdgid_Ws.push_back(ttbarState.W->pdgId());
        tree().pt_Ws.push_back(ttbarState.W->pt());
        tree().eta_Ws.push_back(ttbarState.W->eta());
        tree().phi_Ws.push_back(ttbarState.W->phi());
        tree().mass_Ws.push_back(ttbarState.W->mass());
        tree().energy_Ws.push_back(ttbarState.W->energy());

        if (ttbarState.W_daugther.size() == 2){
          tree().pdgid_dau1.push_back(ttbarState.W_daugther.at(0)->pdgId());
          tree().pt_dau1.push_back(ttbarState.W_daugther.at(0)->pt());
          tree().eta_dau1.push_back(ttbarState.W_daugther.at(0)->eta());
          tree().phi_dau1.push_back(ttbarState.W_daugther.at(0)->phi());
          tree().mass_dau1.push_back(ttbarState.W_daugther.at(0)->mass());
          tree().energy_dau1.push_back(ttbarState.W_daugther.at(0)->energy());

          tree().pdgid_dau2.push_back(ttbarState.W_daugther.at(1)->pdgId());
          tree().pt_dau2.push_back(ttbarState.W_daugther.at(1)->pt());
          tree().eta_dau2.push_back(ttbarState.W_daugther.at(1)->eta());
          tree().phi_dau2.push_back(ttbarState.W_daugther.at(1)->phi());
          tree().mass_dau2.push_back(ttbarState.W_daugther.at(1)->mass());
          tree().energy_dau2.push_back(ttbarState.W_daugther.at(1)->energy());
        }

        tree().pdgid_bjet.push_back(ttbarState.bjet->pdgId());
        tree().pt_bjet.push_back(ttbarState.bjet->pt());
        tree().eta_bjet.push_back(ttbarState.bjet->eta());
        tree().phi_bjet.push_back(ttbarState.bjet->phi());
        tree().mass_bjet.push_back(ttbarState.bjet->mass());
        tree().energy_bjet.push_back(ttbarState.bjet->energy());
      }

     tree.Fill();
}

DEFINE_FWK_MODULE(TTBarGenAnalyzer);
