/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include <iostream>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"

#include "h-tautau/PatProduction/interface/PatVertex.h"

#include "h-tautau/TreeProduction/interface/Muon.h"
#include "h-tautau/TreeProduction/interface/TriggerTools.h"

class MuonBlock : public edm::EDAnalyzer {
public:
    explicit MuonBlock(const edm::ParameterSet& iConfig) :
        _muonInputTag(iConfig.getParameter<edm::InputTag>("muonSrc")),
        _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
        _beamSpotInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
        _beamSpotCorr(iConfig.getParameter<bool>("beamSpotCorr")),
        _muonID(iConfig.getParameter<std::string>("muonID")),
        muonTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { muonTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    edm::InputTag _muonInputTag;
    edm::InputTag _vtxInputTag;
    edm::InputTag _beamSpotInputTag;
    bool _beamSpotCorr;
    std::string _muonID;
    ntuple::MuonTree muonTree;
};

void MuonBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  muonTree.RunId() = iEvent.id().run();
  muonTree.LumiBlock() = iEvent.id().luminosityBlock();
  muonTree.EventId() = iEvent.id().event();

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(_muonInputTag, muons);

  edm::Handle<pat::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByLabel(_beamSpotInputTag, beamSpot);

  if (!muons.isValid()) {
      edm::LogError("MuonBlock") << "Error >> Failed to get pat::Muon Collection for label: "
                                     << _muonInputTag;
      throw std::runtime_error("Failed to get pat::Muon Collection");
  }

  edm::LogInfo("MuonBlock") << "Total # PAT Muons: " << muons->size();
  for (const pat::Muon& patMuon : *muons) {

      // consider only global muons
      if (!patMuon.isGlobalMuon()) continue;

      reco::TrackRef tk  = patMuon.innerTrack();  // tracker segment only
      reco::TrackRef gtk = patMuon.globalTrack();

      muonTree.isTrackerMuon() = patMuon.isTrackerMuon();
      muonTree.isPFMuon() = patMuon.userInt("isPFMuon");

      muonTree.eta()     = patMuon.eta();
      muonTree.phi()     = patMuon.phi();
      muonTree.pt()      = patMuon.pt();
      muonTree.mass() = patMuon.mass();
      muonTree.ptError() = tk->ptError();
      muonTree.charge()  = patMuon.charge();


      double trkd0 = tk->d0();
      double trkdz = tk->dz();
      if (_beamSpotCorr) {
        if (beamSpot.isValid()) {
          trkd0 = -(tk->dxy(beamSpot->position()));
          trkdz = tk->dz(beamSpot->position());
        }
        else{
          edm::LogError("MuonsBlock") << "Error >> Failed to get BeamSpot for label: "
                                      << _beamSpotInputTag;
          throw std::runtime_error("Failed to get BeamSpot for label");
        }
      }
      muonTree.trkD0()  = trkd0;
      muonTree.trkD0Error() = tk->d0Error();
      muonTree.trkDz()  = trkdz;
      muonTree.trkDzError() = tk->dzError();
      muonTree.globalChi2() = patMuon.normChi2();
      muonTree.passID() = patMuon.muonID(_muonID);


      // Vertex association
      double minVtxDist3D = 9999.;
         int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (primaryVertices.isValid()) {
	edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        for (auto vit  = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
          double dxy = tk->dxy(vit->position());
          double dz  = tk->dz(vit->position());
          double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
          if (dist3D < minVtxDist3D) {
            minVtxDist3D = dist3D;
            indexVtx     = int(std::distance(primaryVertices->begin(), vit));
            vertexDistZ  = dz;
          }
        }
      } 
      else {
        edm::LogError("MuonBlock") << "Error >> Failed to get VertexCollection for label: "
                                   << _vtxInputTag;
        throw std::runtime_error("Failed to get VertexCollection for label");
      }
      if(indexVtx >= 0) {
          muonTree.isTightMuon() = patMuon.isTightMuon(primaryVertices->at(indexVtx));
          //          muonTree.isHighPtMuon() = patMuon.isHighPtMuon(primaryVertices->at(indexVtx), reco::improvedTuneP);
          muonTree.isHighPtMuon() = patMuon.isHighPtMuon(primaryVertices->at(indexVtx));
      }
      else {
          muonTree.isTightMuon() = false;
          muonTree.isHighPtMuon() = false;
      }
      muonTree.vtxDist3D() = minVtxDist3D;
      muonTree.vtxIndex()  = indexVtx;
      muonTree.vtxDistZ()  = vertexDistZ;

      // Hit pattern
      const reco::HitPattern& hitp = gtk->hitPattern();  // innerTrack will not provide Muon Hits 
      muonTree.pixHits() = hitp.numberOfValidPixelHits();
      muonTree.trkHits() = hitp.numberOfValidTrackerHits();
      muonTree.muonHits() = hitp.numberOfValidMuonHits();
      muonTree.matches() = patMuon.numberOfMatches();
      muonTree.trackerLayersWithMeasurement() = hitp.trackerLayersWithMeasurement();

      // Isolation
      muonTree.trkIso()  = patMuon.trackIso();
      muonTree.ecalIso() = patMuon.ecalIso();
      muonTree.hcalIso() = patMuon.hcalIso();
      muonTree.hoIso()   = patMuon.isolationR03().hoEt;
      float reliso  = (patMuon.trackIso() + patMuon.ecalIso() + patMuon.hcalIso())/patMuon.pt();
      muonTree.relIso()  = reliso;

      muonTree.pfRelIso() = patMuon.userFloat("PFRelIsoDB04v2");

      // IP information
      muonTree.dB()  = patMuon.dB(pat::Muon::PV2D);
      muonTree.edB() = patMuon.edB(pat::Muon::PV2D);

      muonTree.dB3d()  = patMuon.dB(pat::Muon::PV3D);
      muonTree.edB3d() = patMuon.edB(pat::Muon::PV3D);

      // UW recommendation
      muonTree.isGlobalMuonPromptTight() = muon::isGoodMuon(patMuon, muon::GlobalMuonPromptTight);
      muonTree.isAllArbitrated()         = muon::isGoodMuon(patMuon, muon::AllArbitrated);
      muonTree.nChambers()               = patMuon.numberOfChambers();
      muonTree.nMatches()                = patMuon.numberOfMatches();
      muonTree.nMatchedStations()        = patMuon.numberOfMatchedStations();
      muonTree.stationMask()             = patMuon.stationMask();
      muonTree.stationGapMaskDistance()  = patMuon.stationGapMaskDistance();
      muonTree.stationGapMaskPull()      = patMuon.stationGapMaskPull();
      muonTree.muonID()                  = patMuon.userInt("muonID");

      // Vertex information
      const reco::Candidate::Point& vertex = patMuon.vertex();
      muonTree.vx() = vertex.x();
      muonTree.vy() = vertex.y();
      muonTree.vz() = vertex.z();

      muonTree.idMVA() = patMuon.userFloat("muonIdMVA");
      muonTree.isoRingsMVA() = patMuon.userFloat("muonIsoRingsMVA");
      muonTree.isoRingsRadMVA() = patMuon.userFloat("muonIsoRingsRadMVA");
      muonTree.idIsoCombMVA() = patMuon.userFloat("muonIdIsoCombMVA");

      muonTree.pfRelIso03v1() = patMuon.userFloat("PFRelIso03v1");
      muonTree.pfRelIso03v2() = patMuon.userFloat("PFRelIso03v2");
      muonTree.pfRelIsoDB03v1() = patMuon.userFloat("PFRelIsoDB03v1");
      muonTree.pfRelIsoDB03v2() = patMuon.userFloat("PFRelIsoDB03v2");

      muonTree.pfRelIso04v1() = patMuon.userFloat("PFRelIso04v1");
      muonTree.pfRelIso04v2() = patMuon.userFloat("PFRelIso04v2");
      muonTree.pfRelIsoDB04v1() = patMuon.userFloat("PFRelIsoDB04v1");
      muonTree.pfRelIsoDB04v2() = patMuon.userFloat("PFRelIsoDB04v2");

      muonTree.matchedTriggerPaths() = CollectMatchedTriggerPaths(patMuon);

      muonTree.Fill();
    }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
