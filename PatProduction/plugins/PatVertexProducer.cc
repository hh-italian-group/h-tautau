/*! Definition of PatVertexProducer class which produces analysis-level pat::Vertex objects from reco::Vertexes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "h-tautau/PatProduction/interface/PatVertex.h"

class PatVertexProducer : public edm::EDProducer {
public:
    explicit PatVertexProducer(const edm::ParameterSet&);

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::InputTag _inputTag;
};

PatVertexProducer::PatVertexProducer(const edm::ParameterSet& iConfig)
{
    _inputTag = iConfig.getParameter<edm::InputTag>("inputTag");
    produces<pat::VertexCollection>("");
}

void PatVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByLabel(_inputTag, vertexHandle);
    if (!vertexHandle.isValid()) {
        edm::LogError("PatVertexProducer") << "Error: Failed to get VertexCollection for label: " << _inputTag;
        throw std::runtime_error("Failed to get VertexCollection");
    }

    const reco::VertexCollection& recoVertices = *vertexHandle.product();
    std::auto_ptr<pat::VertexCollection> patVertices(new pat::VertexCollection());
    for(const reco::Vertex& recoVertex : recoVertices) {
        const pat::Vertex patVertex(recoVertex);
        patVertices->push_back(patVertex);
    }

    iEvent.put(patVertices);
}

DEFINE_FWK_MODULE(PatVertexProducer);
