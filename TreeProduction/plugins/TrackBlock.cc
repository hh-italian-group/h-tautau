/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "Math/GenVector/VectorUtil.h"

#include "h-tautau/TreeProduction/interface/Track.h"

class TrackBlock : public edm::EDAnalyzer {
public:
    explicit TrackBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("trackSrc")),
        _beamSpot(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
        trackTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { trackTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;
    edm::InputTag _beamSpot;
    ntuple::TrackTree trackTree;
};

void TrackBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    trackTree.RunId() = iEvent.id().run();
    trackTree.LumiBlock() = iEvent.id().luminosityBlock();
    trackTree.EventId() = iEvent.id().event();

    // read the beam spot
    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByLabel(_beamSpot, beamSpot);

    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByLabel(_inputTag, tracks);

    if (!tracks.isValid()) {
        edm::LogError("TrackBlock") << "Error! Failed to get reco::Track collection, "
                                    << _inputTag;
        throw std::runtime_error("Failed to get reco::Track collection");
    }

    edm::LogInfo("TrackBlock") << "Total # of Tracks: " << tracks->size();
    reco::Track::TrackQuality quality = reco::Track::qualityByName("loose");
    for (const reco::Track& track : *tracks) {
        if (!track.quality(quality)) continue;

        trackTree.eta()         = track.eta();
        trackTree.etaError()    = track.etaError();
        trackTree.theta()       = track.theta();
        trackTree.thetaError()  = track.thetaError();
        trackTree.phi()         = track.phi();
        trackTree.phiError()    = track.phiError();
        trackTree.pt()          = track.pt();
        trackTree.ptError()     = track.ptError();
        trackTree.qoverp()      = track.qoverp();
        trackTree.qoverpError() = track.qoverpError();
        trackTree.charge()      = track.charge();

        trackTree.nValidHits()  = track.numberOfValidHits();
        trackTree.nLostHits()   = track.numberOfLostHits();
        trackTree.validFraction() = track.validFraction();

        const reco::HitPattern& hitp = track.hitPattern();
        trackTree.nValidTrackerHits()  = hitp.numberOfValidTrackerHits();
        trackTree.nValidPixelHits()    = hitp.numberOfValidPixelHits();
        trackTree.nValidStripHits()    = hitp.numberOfValidStripHits();
        trackTree.trackerLayersWithMeasurement() = hitp.trackerLayersWithMeasurement();
        trackTree.pixelLayersWithMeasurement()   = hitp.pixelLayersWithMeasurement();
        trackTree.stripLayersWithMeasurement()   = hitp.stripLayersWithMeasurement();

        trackTree.dxy()       = track.dxy(beamSpot->position());
        trackTree.dxyError()  = track.dxyError();
        trackTree.dz()        = track.dz(beamSpot->position());
        trackTree.dzError()   = track.dzError();

        trackTree.chi2()      = track.chi2();
        trackTree.ndof()      = track.ndof();
        trackTree.vx()        = track.vx();
        trackTree.vy()        = track.vy();
        trackTree.vz()        = track.vz();

        trackTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TrackBlock);
