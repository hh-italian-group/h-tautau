/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <memory>

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "JetMETCorrections/Algorithms/interface/L1FastjetCorrector.h"

class JetUserEmbedder : public edm::EDProducer {
public:
    explicit JetUserEmbedder(const edm::ParameterSet&);
    ~JetUserEmbedder();

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::InputTag jetTag_;
    std::string corrector_;


};

JetUserEmbedder::JetUserEmbedder(const edm::ParameterSet& iConfig)
{

    jetTag_ = iConfig.getParameter<edm::InputTag>("jetTag");
    corrector_ = iConfig.getParameter<std::string>("corrector");


    produces<pat::JetCollection>("");

}


JetUserEmbedder::~JetUserEmbedder()
{

}


void JetUserEmbedder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::JetCollection> jetsHandle;
    iEvent.getByLabel(jetTag_, jetsHandle);
    const pat::JetCollection* jets = jetsHandle.product();


    std::auto_ptr< pat::JetCollection > jetsUserEmbeddedColl( new pat::JetCollection() ) ;

    for(unsigned int i = 0; i < jets->size(); i++) {

        pat::Jet ajet((*jets)[i]);

        const JetCorrector* corrector = JetCorrector::getJetCorrector (corrector_, iSetup);
        const L1FastjetCorrector* l1fastJetCorrector = dynamic_cast<const L1FastjetCorrector*>(corrector);
        const pat::Jet uncorrJet = ajet.correctedJet("Uncorrected");
        float correction = l1fastJetCorrector->correction(uncorrJet,iEvent,iSetup);

        ajet.addUserFloat("correction",correction);

        jetsUserEmbeddedColl->push_back(ajet);

    }

    iEvent.put( jetsUserEmbeddedColl );
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetUserEmbedder);
