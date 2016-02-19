/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "h-tautau/PatProduction/interface/PatVertex.h"

#include "h-tautau/TreeProduction/interface/Vertex.h"

class VertexBlock : public edm::EDAnalyzer {
public:
    explicit VertexBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
        vertexTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { vertexTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;
    ntuple::VertexTree vertexTree;
};

void VertexBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    vertexTree.RunId() = iEvent.id().run();
    vertexTree.LumiBlock() = iEvent.id().luminosityBlock();
    vertexTree.EventId() = iEvent.id().event();

    edm::Handle<pat::VertexCollection> primaryVertices;
    iEvent.getByLabel(_inputTag, primaryVertices);

    if (!primaryVertices.isValid()) {
        edm::LogError("VertexBlock") << "Error >> Failed to get VertexCollection for label: "
                                     << _inputTag;
        throw std::runtime_error("Failed to get VertexCollection");
    }
    edm::LogInfo("VertexBlock") << "Total # Primary Vertices: " << primaryVertices->size();

    for (const pat::Vertex& vertex : *primaryVertices) {
        vertexTree.x() = vertex.x();
        vertexTree.y() = vertex.y();
        vertexTree.z() = vertex.z();
        vertexTree.xErr() = vertex.xError();
        vertexTree.yErr() = vertex.yError();
        vertexTree.zErr() = vertex.zError();
        vertexTree.rho() =  vertex.position().rho();
        vertexTree.chi2() = vertex.chi2();
        vertexTree.ndf() = vertex.ndof();
        vertexTree.ntracks() = vertex.tracksSize();
        vertexTree.ntracksw05() = vertex.nTracks(0.5); // number of tracks in the vertex with weight above 0.5
        vertexTree.isfake() = vertex.isFake();
        vertexTree.isvalid() = vertex.isValid();
        vertexTree.sumPt() = vertex.sumPt();
        vertexTree.sumPtSquared() = vertex.sumPtSquared();

//        edm::LogError("VertexBlock") << vertex.hasRefittedTracks() << " " << vertex.nTracks(0.5) << "\n"
//                                     << vertex.sumPt() << " " << vertex.sumPtSquared() << "\n"
//                                     << vertex.numberOfTracks() << " " << vertex.numberOfWeights();

        vertexTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexBlock);
