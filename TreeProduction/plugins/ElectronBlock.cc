/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>
#include <algorithm>

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "h-tautau/PatProduction/interface/PatVertex.h"
#include "h-tautau/TreeProduction/interface/Electron.h"
#include "h-tautau/TreeProduction/interface/TriggerTools.h"

class ElectronBlock : public edm::EDAnalyzer {
public:
    explicit ElectronBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _bsInputTag(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
        _trkInputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
        _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")),
        _electronInputTag(iConfig.getParameter<edm::InputTag>("electronSrc")),
        electronTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { electronTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _bsInputTag;
    edm::InputTag _trkInputTag;
    edm::InputTag _vtxInputTag;
    edm::InputTag _electronInputTag;
    ntuple::ElectronTree electronTree;
};

void ElectronBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    electronTree.RunId() = iEvent.id().run();
    electronTree.LumiBlock() = iEvent.id().luminosityBlock();
    electronTree.EventId() = iEvent.id().event();

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByLabel(_bsInputTag, beamSpot);

    edm::Handle<pat::VertexCollection> primaryVertices;
    iEvent.getByLabel(_vtxInputTag, primaryVertices);

    edm::Handle<std::vector<pat::Electron> > electrons;
    iEvent.getByLabel(_electronInputTag, electrons);

    if (!electrons.isValid()) {
        edm::LogError("ElectronBlock") << "Error >> Failed to get pat::Electron Collection for label: "
                                       << _electronInputTag;
        throw std::runtime_error("Failed to get pat::Electron Collection");
    }

    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    for (const pat::Electron& patElectron : *electrons) {
      bool hasGsfTrack  = patElectron.gsfTrack().isNonnull() ? true : false;
      reco::GsfTrackRef tk = patElectron.gsfTrack();

      electronTree.ecalDriven() = patElectron.ecalDrivenSeed();
      electronTree.eta() = patElectron.eta();
      electronTree.phi() = patElectron.phi();
      electronTree.pt() = patElectron.pt();
      electronTree.mass() = patElectron.mass();
      electronTree.hasGsfTrack() = hasGsfTrack;
      electronTree.caloEnergy() = patElectron.ecalEnergy();
      electronTree.caloEnergyError() = patElectron.ecalEnergyError();
      electronTree.charge() = patElectron.charge();

      if (hasGsfTrack) {
        electronTree.trackPt() = tk->pt();
        electronTree.trackPtError() = tk->ptError();

    // Hit pattern
    const reco::HitPattern& hitp = tk->hitPattern();
    electronTree.pixHits() = hitp.numberOfValidPixelHits();
    electronTree.trkHits() = hitp.numberOfValidTrackerHits();

        electronTree.nValidHits() = tk->numberOfValidHits();
        electronTree.missingHits() = tk->hitPattern().numberOfHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS);

        electronTree.trkD0() = tk->d0();
        electronTree.trkD0Error() = tk->d0Error();
      }
      // ID variables
      electronTree.hcalOverEcal() = patElectron.hcalOverEcal();
      electronTree.hcalDepth1OverEcal() = patElectron.hcalDepth1OverEcal();
      electronTree.eSuperClusterOverP() = patElectron.eSuperClusterOverP();
      electronTree.sigmaEtaEta() = patElectron.sigmaEtaEta();
      electronTree.sigmaIEtaIEta() = patElectron.sigmaIetaIeta();
      electronTree.deltaPhiTrkSC() = patElectron.deltaPhiSuperClusterTrackAtVtx();
      electronTree.deltaEtaTrkSC() = patElectron.deltaEtaSuperClusterTrackAtVtx();
      electronTree.classification() = patElectron.classification();
      electronTree.e1x5overe5x5() = (patElectron.e5x5() > 0) ? (patElectron.e1x5()/patElectron.e5x5()) : 0;
      electronTree.e2x5overe5x5() = (patElectron.e5x5() > 0) ? (patElectron.e2x5Max()/patElectron.e5x5()) : 0;

      // Iso variables
      electronTree.isoEcal03() = patElectron.dr03EcalRecHitSumEt();
      electronTree.isoHcal03() = patElectron.dr03HcalTowerSumEt();
      electronTree.isoTrk03() = patElectron.dr03TkSumPt();
      electronTree.isoEcal04() = patElectron.dr04EcalRecHitSumEt(); // ecalIso
      electronTree.isoHcal04() = patElectron.dr04HcalTowerSumEt(); // hcalIso
      electronTree.isoTrk04() = patElectron.dr04TkSumPt(); // trackIso
      electronTree.isoRel03() = (patElectron.dr03EcalRecHitSumEt()
                              + patElectron.dr03HcalTowerSumEt()
                              + patElectron.dr03TkSumPt())/patElectron.pt();
      electronTree.isoRel04() = (patElectron.dr04EcalRecHitSumEt()
                              + patElectron.dr04HcalTowerSumEt()
                              + patElectron.dr04TkSumPt())/patElectron.pt();

      // SC associated with electron
      electronTree.scEn() = patElectron.superCluster()->energy();
      electronTree.scEta() = patElectron.superCluster()->eta();
      electronTree.scPhi() = patElectron.superCluster()->phi();
      electronTree.scET() = patElectron.superCluster()->energy()/cosh(patElectron.superCluster()->eta());
      electronTree.scRawEnergy() = patElectron.superCluster()->rawEnergy();

//      electronTree.dist_vec() = patElectron.userFloat("dist");
//      electronTree.dCotTheta() = patElectron.userFloat("dcot");
      electronTree.hasMatchedConversion() = ! patElectron.userInt("antiConv");

      // Vertex association
      double minVtxDist3D = 9999.;
      int indexVtx = -1;
      double vertexDistZ = 9999.;
      if (hasGsfTrack) {
        if (primaryVertices.isValid()) {
      edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
          for (auto vit  = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
            double dxy = tk->dxy(vit->position());
            double dz  = tk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy, 2) + pow(dz, 2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDistZ = dz;
            }
          }
        }
        else {
            edm::LogError("ElectronBlock") << "Error >> Failed to get VertexCollection for label: "
                                           << _vtxInputTag;
            throw std::runtime_error("Failed to get VertexCollection for label");
        }
      }
      // Vertex association variables
      electronTree.vtxDist3D() = minVtxDist3D;
      electronTree.vtxIndex() = indexVtx;
      electronTree.vtxDistZ() = vertexDistZ;

      electronTree.relIso()   = (patElectron.trackIso() + patElectron.ecalIso() + patElectron.hcalIso())/patElectron.pt();

      // PF based isolation
      electronTree.pfRelIso() = patElectron.userFloat("PFRelIsoDB04");

      // PFlow isolation information
      electronTree.chargedHadronIso() = patElectron.chargedHadronIso();
      electronTree.neutralHadronIso() = patElectron.neutralHadronIso();
      electronTree.photonIso() = patElectron.photonIso();

      // IP information
      electronTree.dB() = patElectron.dB(pat::Electron::PV2D);
      electronTree.edB() = patElectron.edB(pat::Electron::PV2D);

      electronTree.dB3d() = patElectron.dB(pat::Electron::PV3D);
      electronTree.edB3d() = patElectron.edB(pat::Electron::PV3D);

      // Bremstrahlung information
      electronTree.nBrems() = patElectron.numberOfBrems();
      electronTree.fbrem() = patElectron.fbrem();

      // MVA
      electronTree.mva() = patElectron.userFloat("mva");
      electronTree.mvaPOGTrig() = patElectron.userFloat("mvaPOGTrig");
      electronTree.mvaPOGNonTrig() = patElectron.userFloat("mvaPOGNonTrig");
      electronTree.mvaPreselection() = (patElectron.userInt("mvaPreselection") ? true : false);
      electronTree.isTriggerElectron() = (patElectron.userInt("isTriggerElectron") ? true : false);

      // MVA Iso
      electronTree.isoMVA() = patElectron.userFloat("eleIsoMVA");

      // Fiducial flag
      unsigned fidFlag = 0;
      if (patElectron.isEB())        fidFlag |= (1 << 0);
      if (patElectron.isEE())        fidFlag |= (1 << 1);
      if (patElectron.isEBEtaGap())  fidFlag |= (1 << 2);
      if (patElectron.isEBPhiGap())  fidFlag |= (1 << 3);
      if (patElectron.isEERingGap()) fidFlag |= (1 << 4);
      if (patElectron.isEEDeeGap())  fidFlag |= (1 << 5);
      if (patElectron.isEBEEGap())   fidFlag |= (1 << 6);
      electronTree.fidFlag() = fidFlag;

      // Vertex information
      const reco::Candidate::Point& vertex = patElectron.vertex();
      electronTree.vx() = vertex.x();
      electronTree.vy() = vertex.y();
      electronTree.vz() = vertex.z();

      electronTree.pfRelIso03v1() = patElectron.userFloat("PFRelIso03v1");
      electronTree.pfRelIso03v2() = patElectron.userFloat("PFRelIso03v2");
      electronTree.pfRelIsoDB03v1() = patElectron.userFloat("PFRelIsoDB03v1");
      electronTree.pfRelIsoDB03v2() = patElectron.userFloat("PFRelIsoDB03v2");
      electronTree.pfRelIsoDB03v3() = patElectron.userFloat("PFRelIsoDB03v3");

      electronTree.pfRelIso04v1() = patElectron.userFloat("PFRelIso04v1");
      electronTree.pfRelIso04v2() = patElectron.userFloat("PFRelIso04v2");
      electronTree.pfRelIsoDB04v1() = patElectron.userFloat("PFRelIsoDB04v1");
      electronTree.pfRelIsoDB04v2() = patElectron.userFloat("PFRelIsoDB04v2");
      electronTree.pfRelIsoDB04v3() = patElectron.userFloat("PFRelIsoDB04v3");

      // 2012
      electronTree.pfRelIso03() = patElectron.userFloat("PFRelIso03");
      electronTree.pfRelIso04() = patElectron.userFloat("PFRelIso04");
      electronTree.pfRelIsoDB03() = patElectron.userFloat("PFRelIsoDB03");
      electronTree.pfRelIsoDB04() = patElectron.userFloat("PFRelIsoDB04");

      electronTree.matchedTriggerPaths() = CollectMatchedTriggerPaths(patElectron);

      electronTree.Fill();
    }
}

DEFINE_FWK_MODULE(ElectronBlock);
