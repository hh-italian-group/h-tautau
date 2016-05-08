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
//
// class declaration
//

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

      edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
      edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
      edm::EDGetTokenT<edm::View<reco::GenJet> > genJetToken_;
};

TTBarGenAnalyzer::TTBarGenAnalyzer(const edm::ParameterSet& iConfig):
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
  genJetToken_(consumes< edm::View<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("genJet")))   
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

        for(size_t i=0; i<pruned->size();i++){
                // if(abs((*pruned)[i].pdgId()) > 500 && abs((*pruned)[i].pdgId()) <600)
                if(abs((*pruned)[i].pdgId()) == 6 ){
                        const Candidate * top = &(*pruned)[i];
                        if ((*pruned)[i].status() == 62){
                        std::cout << "PdgID: " << top->pdgId() << " pt " << top->pt() << " eta: " << top->eta() << " phi: " << top->phi() << std::endl;
                        std::cout << " Number of Daughters:  " << top->numberOfDaughters() << "  Status: " << top->status() << std::endl;
                        for (size_t i=0; i<top->numberOfDaughters(); ++i){
                           const Candidate * daughter = top->daughter(i);
                           std::cout << "\t\t PdgID: " << daughter->pdgId() << " pt " << daughter->pt() << " eta: " << daughter->eta() << " phi: " << daughter->phi() << std::endl;
                           
                           if ( std::abs(daughter->pdgId()) == 24){
                              auto W = findDecayParticle(daughter);
                              for (size_t i=0; i<W->numberOfDaughters(); ++i){
                                 const Candidate * Wdaughter = W->daughter(i);
                                 std::cout << "\t\t\t PdgID: " << Wdaughter->pdgId() << " pt " << Wdaughter->pt() << " eta: " << Wdaughter->eta() << " phi: " << Wdaughter->phi() << std::endl;
                              }
                           }
                           if ( std::abs(daughter->pdgId()) == 5 ){                               
                              for (size_t i=0; i<genJet->size();i++){
                                 const GenJet * jet = &(*genJet)[i];
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
                        
                       }
                }

        }


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
}

DEFINE_FWK_MODULE(TTBarGenAnalyzer);