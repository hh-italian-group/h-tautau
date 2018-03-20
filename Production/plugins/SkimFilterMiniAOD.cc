/*! Definition of the class that applies skimming of the miniAOD files.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

class SkimFilterMiniAOD : public edm::EDFilter {
public:
    explicit SkimFilterMiniAOD(const edm::ParameterSet& iConfig)
        : vertexTag(mayConsume<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
          muonTag(mayConsume<pat::MuonCollection >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
          electronTag(mayConsume<pat::ElectronCollection>(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
          tauTag(mayConsume<pat::TauCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))) {}

private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);

private:
//    const edm::EDGetTokenT<reco::VertexCollection> vertexTag;
    edm::EDGetToken vertexTag;
    edm::EDGetToken muonTag;
    edm::EDGetToken electronTag;
    edm::EDGetToken tauTag;
};

bool SkimFilterMiniAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::VertexCollection> hVertices;
    iEvent.getByToken(vertexTag, hVertices);
    const reco::VertexCollection& vertices = *hVertices;

    edm::Handle<pat::MuonCollection> hMuons;
    iEvent.getByToken(muonTag, hMuons);
    const pat::MuonCollection& muons = *hMuons.product();

    edm::Handle<pat::ElectronCollection> hElectrons;
    iEvent.getByToken(electronTag, hElectrons);
    const pat::ElectronCollection& electrons = *hElectrons.product();

    edm::Handle<pat::TauCollection> hTaus;
    iEvent.getByToken(tauTag, hTaus);
    const pat::TauCollection& taus = *hTaus.product();

    bool haveGoodVertex = false;
    unsigned nGoodTaus = 0;

  //  auto vtx = vertices.at(0);

    for (const reco::Vertex& vertex : vertices){
        if(vertex.ndof() > 4 && std::abs(vertex.z()) < 24 && vertex.position().rho() < 2){
            haveGoodVertex = true;
            break;
        }
    }
    if (!haveGoodVertex) {
        return false;}
  

    for(const pat::Tau& tau : taus) {
        if(tau.pt() > 18 && std::abs(tau.eta()) < 2.3  /*&& tau.tauID("decayModeFinding") > 0.5*/)
            ++nGoodTaus;
        if(nGoodTaus >= 2) return true;
    }
    if(!nGoodTaus) return false;
  

    for(const pat::Electron& electron : electrons) {
        if(electron.pt() > 20.0 && std::abs(electron.eta()) < 2.0) {
            std::cout<<" Ele Missing Hits =  "<<
                       electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS)<<std::endl;
            return true;}
    }

    for(const pat::Muon& muon : muons) {
//        bool muonIP = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.045 &&
//                      fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.2;
        if(muon.pt() > 16.0 && std::abs(muon.eta())<2.4 && muon.isLooseMuon()) {
//            float iso = (muon.pfIsolationR03().sumChargedHadronPt + std::max(
//                       muon.pfIsolationR03().sumNeutralHadronEt +
//                       muon.pfIsolationR03().sumPhotonEt -
//                       0.5 * muon.pfIsolationR03().sumPUPt, 0.0)) / muon.pt();

//            std::cout<<" Muon Isolament =  "<<iso<<std::endl;

            return true;
        }
    }

    return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimFilterMiniAOD);
