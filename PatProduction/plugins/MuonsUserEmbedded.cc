/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <memory>

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

namespace {
template<typename T>
bool isValidRef(const edm::Ref<T>& ref)
{
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
}
}

class MuonsUserEmbedded : public edm::EDProducer {
public:
    explicit MuonsUserEmbedded(const edm::ParameterSet&);

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::InputTag muonTag_;
    edm::InputTag vertexTag_;

    reco::isodeposit::AbsVetos vetos2010Charged_;
    reco::isodeposit::AbsVetos vetos2010Neutral_;
    reco::isodeposit::AbsVetos vetos2010Photons_;

    reco::isodeposit::AbsVetos vetos2011Charged_;
    reco::isodeposit::AbsVetos vetos2011Neutral_;
    reco::isodeposit::AbsVetos vetos2011Photons_;
};


MuonsUserEmbedded::MuonsUserEmbedded(const edm::ParameterSet& iConfig)
{
    muonTag_ = iConfig.getParameter<edm::InputTag>("muonTag");
    vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");

    produces<pat::MuonCollection>("");
}


void MuonsUserEmbedded::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::MuonCollection> muonsHandle;
    iEvent.getByLabel(muonTag_, muonsHandle);
    const pat::MuonCollection* muons = muonsHandle.product();

//  edm::Handle<reco::MuonCollection> recoMuonsHandle;
//  iEvent.getByLabel("muons",recoMuonsHandle);
    //const reco::MuonCollection* recoMuons = recoMuonsHandle.product();

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByLabel(vertexTag_, vertexHandle);
    const reco::VertexCollection* vertexes = vertexHandle.product();

    edm::Handle<reco::PFCandidateCollection> pfHandle;
    iEvent.getByLabel("particleFlow", pfHandle);
    if( !pfHandle.isValid() )
        edm::LogError("DataNotAvailable")
                << "No pf particles label available \n";
    const reco::PFCandidateCollection* pfCandidates = pfHandle.product();

    std::auto_ptr< pat::MuonCollection > muonsUserEmbeddedColl( new pat::MuonCollection() ) ;

    for(unsigned int i = 0; i < muons->size(); i++) {
        pat::Muon aMuon( (*muons)[i] );

        //const reco::Muon* aRecoMuon = 0;
        /*    for(unsigned int j = 0; j < recoMuons->size(); j++){

               if( Geom::deltaR( (*recoMuons)[j].p4() , aMuon.p4()) < 1e-03 ) {
        aRecoMuon = &((*recoMuons)[j]);
        std::cout << "Match to recoMuon" <<aRecoMuon->pt()<< std::endl;
        }
        }*/

        //std::cout << "@@@@@@@@ recoMuon: " << aMuon.px() << ", " << aMuon.py() << ", " << aMuon.pz() << std::endl;
        int isPFMuon = 0;
        for(unsigned int j = 0; j < pfCandidates->size(); j++) {
            if( (*pfCandidates)[j].particleId() == reco::PFCandidate::mu ) {
                reco::MuonRef muonRefToPFMuon = (*pfCandidates)[j].muonRef();
                //if( muonRefToPFMuon.isNonnull() ){
                //std::cout << j << ": muonRefToPFMuon: " << muonRefToPFMuon->px() << ", "
                //    << muonRefToPFMuon->py() << ", " << muonRefToPFMuon->pz() << std::endl;
                //}
                if( muonRefToPFMuon.isNonnull() &&
                        Geom::deltaR( muonRefToPFMuon->p4() , aMuon.p4()) < 1e-04 &&
                        (muonRefToPFMuon->isGlobalMuon() || muonRefToPFMuon->isTrackerMuon() ) )
                    isPFMuon = 1;
            }
        }

        double dxyWrtPV =  -99.;
        double dzWrtPV  =  -99.;

        if(vertexes->size() != 0 && aMuon.isGlobalMuon()) {
            dxyWrtPV = (aMuon.globalTrack())->dxy( (*vertexes)[0].position() ) ;
            dzWrtPV  = (aMuon.globalTrack())->dz( (*vertexes)[0].position() ) ;
        } else if (vertexes->size() != 0 && aMuon.isTrackerMuon()) {
            dxyWrtPV = (aMuon.innerTrack())->dxy( (*vertexes)[0].position() ) ;
            dzWrtPV  = (aMuon.innerTrack())->dz( (*vertexes)[0].position() ) ;
        }

        aMuon.addUserFloat("dxyWrtPV", dxyWrtPV);
        aMuon.addUserFloat("dzWrtPV", dzWrtPV);
        aMuon.addUserInt("isPFMuon", isPFMuon);

        int isGlobalMuon_ = 0;
        int isTrackerMuon_ = 0;
        int globalMuonPromptTight_ = 0;
        int allArbitrated_ = 0;
        double normalizeChi2_ = 99;
        double ptError_ = 99;
        int innerTrackHitPattern_ = -99;
        int numberOfValidPixelHits_ = -99;
        int numMuonStations_ = -99;
        int numberOfMatches_ = -99;

        isGlobalMuon_ = aMuon.isGlobalMuon();
        isTrackerMuon_ = aMuon.isTrackerMuon();

        bool validRefGlob = isValidRef(aMuon.globalTrack());
        bool validRefInn = isValidRef(aMuon.innerTrack());

        if(validRefGlob && validRefInn) {

            globalMuonPromptTight_ = muon::isGoodMuon(aMuon, muon::GlobalMuonPromptTight);
            allArbitrated_ =  muon::isGoodMuon(aMuon, muon::AllArbitrated);

            normalizeChi2_ = aMuon.globalTrack()->normalizedChi2();
            ptError_ = (aMuon.innerTrack()->ptError()) / (aMuon.innerTrack()->pt());

            const reco::HitPattern& innerTrackHitPattern = aMuon.innerTrack()->hitPattern();
            innerTrackHitPattern_ = innerTrackHitPattern.numberOfValidTrackerHits();
            numberOfValidPixelHits_ = innerTrackHitPattern.numberOfValidPixelHits();

            int numMuonStations = 0;
            unsigned int theStationMask = (unsigned int)aMuon.stationMask(reco::Muon::SegmentAndTrackArbitration);
            for( int i = 0; i < 8; ++i ) { // eight stations, eight bits

                if ( theStationMask & (1 << i) ) ++numMuonStations;

            }

            numMuonStations_ = numMuonStations;
            numberOfMatches_ = aMuon.numberOfMatches();

        }

        int muonID = isGlobalMuon_ && isTrackerMuon_ && globalMuonPromptTight_ && allArbitrated_
                && (fabs(dxyWrtPV) < 0.02) && (fabs(dzWrtPV) < 0.2) && (normalizeChi2_ < 10) && (ptError_ < 0.1)
                && (innerTrackHitPattern_ >= 10) && (numberOfValidPixelHits_ >= 1) && (numMuonStations_ >= 2)
                && (numberOfMatches_ >= 1);

        //std::cout <<"isGlobal: "<< isGlobalMuon_<<" isTracker: "<<isTrackerMuon_
        //          <<" globalMuonPromptTight: "<<globalMuonPromptTight_<<" allArbitrated"<< allArbitrated_
        //          <<" dxyWrtPV: "<<dxyWrtPV<<" dzWrtPV: "<<dzWrtPV<<" normalizeChi2: "<<normalizeChi2_
        //          <<" ptError: "<<ptError_<<" innerTrackHitPattern: "<<innerTrackHitPattern_
        //          <<" numberOfValidPixelHits: "<<numberOfValidPixelHits_<<" numMuonStations: "<<numMuonStations_
        //          <<" numberOfMatches: "<<numberOfMatches_<<std::endl;

        //std::cout<<"Total ID: "<<muonID<<std::endl;

        aMuon.addUserInt("isGlobalMuon", isGlobalMuon_);
        aMuon.addUserInt("isTrackerMuon", isTrackerMuon_);
        aMuon.addUserInt("globalMuonPromptTight", globalMuonPromptTight_);
        aMuon.addUserInt("allArbitrated", allArbitrated_);
        aMuon.addUserFloat("normalizeChi2", normalizeChi2_);
        aMuon.addUserFloat("ptError", ptError_);
        aMuon.addUserInt("innerTrackHitPattern", innerTrackHitPattern_);
        aMuon.addUserInt("numberOfValidPixelHits", numberOfValidPixelHits_);
        aMuon.addUserInt("numMuonStations", numMuonStations_);
        aMuon.addUserInt("numberOfMatches", numberOfMatches_);
        aMuon.addUserInt("muonID", muonID);

        // iso deposits
        reco::isodeposit::AbsVetos vetos2010Charged;
        reco::isodeposit::AbsVetos vetos2010Neutral;
        reco::isodeposit::AbsVetos vetos2010Photons;
        reco::isodeposit::AbsVetos vetos2011Charged;
        reco::isodeposit::AbsVetos vetos2011Neutral;
        reco::isodeposit::AbsVetos vetos2011Photons;

        vetos2010Charged.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aMuon.eta(), aMuon.phi()), 0.01));
        vetos2010Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));
        vetos2010Neutral.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aMuon.eta(), aMuon.phi()), 0.08));
        vetos2010Neutral.push_back(new reco::isodeposit::ThresholdVeto(1.0));
        vetos2010Photons.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aMuon.eta(), aMuon.phi()), 0.05));
        vetos2010Photons.push_back(new reco::isodeposit::ThresholdVeto(1.0));

        vetos2011Charged.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aMuon.eta(), aMuon.phi()), 0.0001));
        vetos2011Charged.push_back(new reco::isodeposit::ThresholdVeto(0.0));
        vetos2011Neutral.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aMuon.eta(), aMuon.phi()), 0.01));
        vetos2011Neutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));
        vetos2011Photons.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aMuon.eta(), aMuon.phi()), 0.01));
        vetos2011Photons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

        float chIso03v1 =
            aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2010Charged).first;
        float nhIso03v1 =
            aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2010Neutral).first;
        float phIso03v1 =
            aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2010Photons).first;
        float nhIsoPU03v1 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2010Neutral).first;
        float phIsoPU03v1 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2010Photons).first;

        float chIso04v1 =
            aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2010Charged).first;
        float nhIso04v1 =
            aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2010Neutral).first;
        float phIso04v1 =
            aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2010Photons).first;
        float nhIsoPU04v1 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2010Neutral).first;
        float phIsoPU04v1 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2010Photons).first;

        float chIso03v2 =
            aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2011Charged).first;
        float nhIso03v2 =
            aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2011Neutral).first;
        float phIso03v2 =
            aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2011Photons).first;
        float nhIsoPU03v2 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2011Neutral).first;
        float phIsoPU03v2 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2011Photons).first;

        // User1Iso corresponds in the python (see patMuons.py) to muPFIsoDepositChargedAllPFIso -> all noPU CH+MU+E
        float allChIso04v2 =
            aMuon.isoDeposit(pat::User1Iso)->depositAndCountWithin(0.4, vetos2011Charged).first;
        float chIso04v2 =
            aMuon.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2011Charged).first;
        float nhIso04v2 =
            aMuon.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2011Neutral).first;
        float phIso04v2 =
            aMuon.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2011Photons).first;
        float nhIsoPU04v2 =
            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011Neutral).first;
//        float phIsoPU04v2 =
//            aMuon.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011Photons).first;


        aMuon.addUserFloat("PFRelIso04v1", (chIso04v1 + nhIso04v1 + phIso04v1) / aMuon.pt());
        aMuon.addUserFloat("PFRelIso03v1", (chIso03v1 + nhIso03v1 + phIso03v1) / aMuon.pt());
        aMuon.addUserFloat("PFRelIsoDB04v1",
            (chIso04v1 + std::max(nhIso04v1 + phIso04v1 - 0.5 * 0.5 * (nhIsoPU04v1 + phIsoPU04v1), 0.0)) / aMuon.pt());
        aMuon.addUserFloat("PFRelIsoDB03v1",
            (chIso03v1 + std::max(nhIso03v1 + phIso03v1 - 0.5 * 0.5 * (nhIsoPU03v1 + phIsoPU03v1), 0.0)) / aMuon.pt());

        aMuon.addUserFloat("PFRelIso04v2", (chIso04v2 + nhIso04v2 + phIso04v2) / aMuon.pt());
        aMuon.addUserFloat("PFRelIso03v2", (chIso03v2 + nhIso03v2 + phIso03v2) / aMuon.pt());
        // new definition - Muon POG
//        aMuon.addUserFloat("PFRelIsoDB04v2",
//            (chIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * 0.5 * (nhIsoPU04v2 + phIsoPU04v2), 0.0)) / aMuon.pt());
        //legacy definition
        aMuon.addUserFloat("PFRelIsoDB04v2",
            (allChIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * (nhIsoPU04v2), 0.0)) / aMuon.pt());
//        if (iEvent.id().event() == 11370491 && iEvent.run() == 190706 && iEvent.luminosityBlock() == 15){
//            std::cout << "event: " << iEvent.id() << ", allChIso04v2: " << allChIso04v2 << ", chIso04v2: " << chIso04v2
//                      << ", nhIso04v2: " << nhIso04v2 << ", phIso04v2: " << phIso04v2 << ", nhIsoPU04v2: " << nhIsoPU04v2
//                      << ", phIsoPU04v2: " << phIsoPU04v2 << ", sum: " <<
//                         (allChIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * (nhIsoPU04v2), 0.0)) <<
//                         ", pt Mu: " << aMuon.pt() << ", eta Mu: " << aMuon.eta() << ", phi Mu: " << aMuon.phi() <<
//                         ", iso: " << (allChIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * (nhIsoPU04v2), 0.0)) / aMuon.pt() <<
//                         std::endl;
//        }
        aMuon.addUserFloat("PFRelIsoDB03v2",
            (chIso03v2 + std::max(nhIso03v2 + phIso03v2 - 0.5 * 0.5 * (nhIsoPU03v2 + phIsoPU03v2), 0.0)) / aMuon.pt());

        // cleaning
        for(unsigned int i = 0; i < vetos2010Charged.size(); i++) {
            delete vetos2010Charged[i];
        }
        for(unsigned int i = 0; i < vetos2010Neutral.size(); i++) {
            delete vetos2010Neutral[i];
            delete vetos2010Photons[i];
        }
        for(unsigned int i = 0; i < vetos2011Charged.size(); i++) {
            delete vetos2011Charged[i];
        }
        for(unsigned int i = 0; i < vetos2011Neutral.size(); i++) {
            delete vetos2011Neutral[i];
            delete vetos2011Photons[i];
        }

        aMuon.addUserFloat("isInRun", iEvent.run());

        muonsUserEmbeddedColl->push_back(aMuon);

    }

    iEvent.put( muonsUserEmbeddedColl );
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonsUserEmbedded);
