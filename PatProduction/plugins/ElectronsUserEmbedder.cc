/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <memory>

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimator.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/plugins/TransientTrackBuilderESProducer.h"

class ElectronsUserEmbedder : public edm::EDProducer {
public:
    explicit ElectronsUserEmbedder(const edm::ParameterSet&);
    ~ElectronsUserEmbedder();

private:
    virtual void produce(edm::Event&, const edm::EventSetup&);

    edm::InputTag electronTag_;
    edm::InputTag vertexTag_;
    bool isMC_;
    bool doMVAPOG_;

    EGammaMvaEleEstimator* myMVATrig_;
    EGammaMvaEleEstimator* myMVANonTrig_;

};

ElectronsUserEmbedder::ElectronsUserEmbedder(const edm::ParameterSet& iConfig)
{

    electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
    vertexTag_   = iConfig.getParameter<edm::InputTag>("vertexTag");
    isMC_        = iConfig.getParameter<bool>("isMC");
    doMVAPOG_    = iConfig.getParameter<bool>("doMVAPOG");

    edm::FileInPath inputFileName0v2 = iConfig.getParameter<edm::FileInPath>("inputFileName0v2");
    edm::FileInPath inputFileName1v2 = iConfig.getParameter<edm::FileInPath>("inputFileName1v2");
    edm::FileInPath inputFileName2v2 = iConfig.getParameter<edm::FileInPath>("inputFileName2v2");
    edm::FileInPath inputFileName3v2 = iConfig.getParameter<edm::FileInPath>("inputFileName3v2");
    edm::FileInPath inputFileName4v2 = iConfig.getParameter<edm::FileInPath>("inputFileName4v2");
    edm::FileInPath inputFileName5v2 = iConfig.getParameter<edm::FileInPath>("inputFileName5v2");

    edm::FileInPath inputFileName0v3 = iConfig.getParameter<edm::FileInPath>("inputFileName0v3");
    edm::FileInPath inputFileName1v3 = iConfig.getParameter<edm::FileInPath>("inputFileName1v3");
    edm::FileInPath inputFileName2v3 = iConfig.getParameter<edm::FileInPath>("inputFileName2v3");
    edm::FileInPath inputFileName3v3 = iConfig.getParameter<edm::FileInPath>("inputFileName3v3");
    edm::FileInPath inputFileName4v3 = iConfig.getParameter<edm::FileInPath>("inputFileName4v3");
    edm::FileInPath inputFileName5v3 = iConfig.getParameter<edm::FileInPath>("inputFileName5v3");

    if(doMVAPOG_) {

        Bool_t manualCat = true;

        std::vector<string> myManualCatWeigthsTrig;
        myManualCatWeigthsTrig.push_back(inputFileName0v2.fullPath().data());
        myManualCatWeigthsTrig.push_back(inputFileName1v2.fullPath().data());
        myManualCatWeigthsTrig.push_back(inputFileName2v2.fullPath().data());
        myManualCatWeigthsTrig.push_back(inputFileName3v2.fullPath().data());
        myManualCatWeigthsTrig.push_back(inputFileName4v2.fullPath().data());
        myManualCatWeigthsTrig.push_back(inputFileName5v2.fullPath().data());

        myMVATrig_ = new EGammaMvaEleEstimator();
        myMVATrig_->initialize("BDT",
                               EGammaMvaEleEstimator::kTrig,
                               manualCat,
                               myManualCatWeigthsTrig);

        std::vector<string> myManualCatWeigthsNonTrig;
        myManualCatWeigthsNonTrig.push_back(inputFileName0v3.fullPath().data());
        myManualCatWeigthsNonTrig.push_back(inputFileName1v3.fullPath().data());
        myManualCatWeigthsNonTrig.push_back(inputFileName2v3.fullPath().data());
        myManualCatWeigthsNonTrig.push_back(inputFileName3v3.fullPath().data());
        myManualCatWeigthsNonTrig.push_back(inputFileName4v3.fullPath().data());
        myManualCatWeigthsNonTrig.push_back(inputFileName5v3.fullPath().data());

        myMVANonTrig_ = new EGammaMvaEleEstimator();
        myMVANonTrig_->initialize("BDT",
                                  EGammaMvaEleEstimator::kNonTrig,
                                  manualCat,
                                  myManualCatWeigthsNonTrig);

    }

    produces<pat::ElectronCollection>("");

}


ElectronsUserEmbedder::~ElectronsUserEmbedder()
{
    if(doMVAPOG_) {
        delete myMVATrig_; //delete fMVADaniele_;
        delete myMVANonTrig_;
    }
}


void ElectronsUserEmbedder::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<pat::ElectronCollection> electronsHandle;
    iEvent.getByLabel(electronTag_, electronsHandle);
    const pat::ElectronCollection* electrons = electronsHandle.product();

    edm::Handle<reco::VertexCollection> vertexHandle;
    iEvent.getByLabel(vertexTag_, vertexHandle);
    const reco::VertexCollection* vertexes = vertexHandle.product();

    edm::Handle<reco::BeamSpot> bsHandle;
    iEvent.getByLabel("offlineBeamSpot", bsHandle);
    const reco::BeamSpot &thebs = *bsHandle.product();

    edm::Handle<reco::ConversionCollection> hConversions;
    iEvent.getByLabel("allConversions", hConversions);

    edm::Handle<DcsStatusCollection> dcsHandle;
    iEvent.getByLabel("scalersRawToDigi", dcsHandle);
//    float evt_bField;

    if (!isMC_) {
        // scale factor = 3.801/18166.0 which are
        // average values taken over a stable two
        // week period
//        float currentToBFieldScaleFactor = 2.09237036221512717e-04;
//        float current = -9999 / currentToBFieldScaleFactor;
//        if( dcsHandle.isValid() && (*dcsHandle).size() > 0 ) {
//            current = (*dcsHandle)[0].magnetCurrent();
//        }

//        evt_bField = current * currentToBFieldScaleFactor;
//    } else {
//        edm::ESHandle<MagneticField> magneticField;
//        iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
//        evt_bField = magneticField->inTesla(GlobalPoint(0., 0., 0.)).z();
    }

    // PFCandidateMap
    //edm::Handle<edm::ValueMap<reco::PFCandidatePtr> >ValMapH;
    //iEvent.getByLabel(edm::InputTag("particleFlow","electrons"),ValMapH); // The InputTag is 'particleFlow:electrons'
    //const edm::ValueMap<reco::PFCandidatePtr> & myValMap(*ValMapH);

    //Get the CTF tracks
    edm::Handle<reco::TrackCollection> tracks_h;
    iEvent.getByLabel("generalTracks", tracks_h);

    //get GSF Tracks
    edm::Handle<reco::GsfTrackCollection> gsftracks_h;
    iEvent.getByLabel("electronGsfTracks", gsftracks_h);

    edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", hTransientTrackBuilder);
//    const TransientTrackBuilder *transientTrackBuilder = hTransientTrackBuilder.product();

    std::auto_ptr< pat::ElectronCollection > electronsUserEmbeddedColl( new pat::ElectronCollection() ) ;

    for(unsigned int i = 0; i < electrons->size(); i++) {

        pat::Electron aElectron( (*electrons)[i] );
        const reco::GsfElectron* aGsf = static_cast<reco::GsfElectron*>(&aElectron);

        //unsigned int position = 0;
        //bool matchedToGsf = false;
        //for(unsigned int j = 0; j < gsfElectrons->size(); j++){
        //if( Geom::deltaR((*gsfElectrons)[j].p4(), aGsf->p4()) < 1e-03 ){
        //matchedToGsf = true;
        //position = j;
        //}
        //}
        //reco::GsfElectronRef aGsfElectronRef(gsfElectronsHandle, position);
        //const reco::PFCandidatePtr pfElePtr(myValMap[aGsfElectronRef]);
        int pfId = 1; //pfElePtr.isNonnull();

        const reco::Track *el_track = (const reco::Track*)((aElectron).gsfTrack().get());
        const reco::HitPattern& p_inner = el_track->hitPattern();
        float nHits = p_inner.numberOfHits(reco::HitPattern::HitCategory::TRACK_HITS);

        ConversionFinder convFinder;
//        vector<ConversionInfo> v_convInfos = convFinder.getConversionInfos(*(aElectron.core()), tracks_h,
//                                                                           gsftracks_h, evt_bField);
//        std::cout << "before conversion" << std::endl;
//        ConversionInfo convInfo  = convFinder.getConversionInfo(*aGsf, tracks_h, gsftracks_h, evt_bField);
//        std::cout << "after conversion" << std::endl;
//        double els_conv_dist     = convInfo.dist();
//        double els_conv_dcot     = convInfo.dcot();
        //double els_conv_radius = convInfo.radiusOfConversion();
        //math::XYZPoint els_conv_Point = convInfo.pointOfConversion();
        //    std::cout<<els_conv_Point<<std::endl;
//        reco::TrackRef els_conv_ctfRef = convInfo.conversionPartnerCtfTk();
//        reco::GsfTrackRef els_conv_gsfRef = convInfo.conversionPartnerGsfTk();
        //double els_conv_delMissHits = convInfo.deltaMissingHits();

        float dPhi  = aElectron.deltaPhiSuperClusterTrackAtVtx();
        float dEta  = aElectron.deltaEtaSuperClusterTrackAtVtx();
        float sihih = aElectron.sigmaIetaIeta();
        float HoE   = aElectron.hadronicOverEm();
        //cout << "dEta " << dEta << " dPhi " << dPhi << " -- dcot " << els_conv_dcot << " -- nHits " << nHits << endl;

        aElectron.addUserFloat("nHits", nHits);
//        aElectron.addUserFloat("dist", fabs(els_conv_dist));
//        aElectron.addUserFloat("dcot", fabs(els_conv_dcot));
        aElectron.addUserFloat("dPhi", fabs(dPhi));
        aElectron.addUserFloat("dEta", fabs(dEta));
        aElectron.addUserFloat("sihih", sihih);
        aElectron.addUserFloat("HoE", HoE);

        int passconversionveto =
            int(!ConversionTools::hasMatchedConversion(*aGsf, hConversions, thebs.position(), true, 2.0, 1e-06, 0));
        aElectron.addUserInt("antiConv", passconversionveto);

        //edm::Handle< EcalRecHitCollection > reducedEBRecHits;
        //edm::Handle< EcalRecHitCollection > reducedEERecHits;
        //iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEB"), reducedEBRecHits );
        //iEvent.getByLabel( edm::InputTag("reducedEcalRecHitsEE"), reducedEERecHits ) ;
        //const EcalRecHitCollection * reducedRecHits = 0 ;
        //if (aElectron.isEB())
        //  reducedRecHits = reducedEBRecHits.product() ;
        //else
        //  reducedRecHits = reducedEERecHits.product() ;

        double dxyWrtPV = -99.;
        double dzWrtPV =  -99.;

        if(vertexes->size() != 0 && (aElectron.gsfTrack()).isNonnull() ) {
            dxyWrtPV = (aElectron.gsfTrack())->dxy( (*vertexes)[0].position() ) ;
            dzWrtPV  = (aElectron.gsfTrack())->dz(  (*vertexes)[0].position() ) ;
        } else if (vertexes->size() != 0 && (aElectron.track()).isNonnull() ) {
            dxyWrtPV = (aElectron.track())->dxy( (*vertexes)[0].position() ) ;
            dzWrtPV  = (aElectron.track())->dz(  (*vertexes)[0].position() ) ;
        }

        aElectron.addUserFloat("dxyWrtPV", dxyWrtPV);
        aElectron.addUserFloat("dzWrtPV", dzWrtPV);

//        EcalClusterLazyTools lazyTools(iEvent, iSetup,
//                                       edm::InputTag("reducedEcalRecHitsEB"),
//                                       edm::InputTag("reducedEcalRecHitsEE"));

        int myTrigPresel = 0;
        if(fabs(aGsf->superCluster()->eta()) < 1.485) {
            if(aGsf->sigmaIetaIeta() < 0.014 &&
                    aGsf->hadronicOverEm() < 0.15 &&
                    aGsf->dr03TkSumPt() / aGsf->pt() < 0.2 &&
                    aGsf->dr03EcalRecHitSumEt() / aGsf->pt() < 0.2 &&
                    aGsf->dr03HcalTowerSumEt() / aGsf->pt() < 0.2 &&
                    aGsf->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS) == 0)
                myTrigPresel = 1;
        } else {
            if(aGsf->sigmaIetaIeta() < 0.035 &&
                    aGsf->hadronicOverEm() < 0.10 &&
                    aGsf->dr03TkSumPt() / aGsf->pt() < 0.2 &&
                    aGsf->dr03EcalRecHitSumEt() / aGsf->pt() < 0.2 &&
                    aGsf->dr03HcalTowerSumEt() / aGsf->pt() < 0.2 &&
                    aGsf->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::HitCategory::MISSING_INNER_HITS) == 0)
                myTrigPresel = 1;
        }

        float mva2 = -99;
        float mva3 = -99;
        int mvaPreselection = passconversionveto && nHits <= 0 && dxyWrtPV < 0.02 && dzWrtPV < 0.1 &&
                              ((aElectron.isEB() && sihih < 0.01
                                && fabs(dEta) < 0.007
                                && fabs(dPhi) < 0.15
                                && HoE < 0.12
                                && aElectron.dr03TkSumPt() / aElectron.pt() < 0.20
                                && (TMath::Max(aElectron.dr03EcalRecHitSumEt() - 1.0, 0.0)) / aElectron.pt() < 0.20
                                && aElectron.dr03HcalTowerSumEt() / aElectron.pt() < 0.20

                               ) ||
                               (aElectron.isEE() && sihih < 0.03
                                && fabs(dEta) < 0.009
                                && fabs(dPhi) < 0.10
                                && HoE < 0.10
                                && aElectron.dr03TkSumPt() / aElectron.pt() < 0.20
                                && (TMath::Max(aElectron.dr03EcalRecHitSumEt() - 1.0, 0.0)) / aElectron.pt() < 0.20
                                && aElectron.dr03HcalTowerSumEt() / aElectron.pt() < 0.20
                               ));

        if(doMVAPOG_) {
            //mva2 = fMVADaniele_->mva(*aGsf, vertexes->size());
//            mva2 = myMVATrig_->mvaValue( *aGsf , (*vertexes)[0], *transientTrackBuilder, lazyTools, false);
//            mva3 = myMVANonTrig_->mvaValue( *aGsf , (*vertexes)[0], *transientTrackBuilder, lazyTools, false);
            //cout << mva2 << endl;
        }

        aElectron.addUserInt("mvaPreselection", mvaPreselection);
        aElectron.addUserInt("isTriggerElectron", myTrigPresel);
        aElectron.addUserFloat("mvaPOGTrig"   , mva2);
        aElectron.addUserFloat("mvaPOGNonTrig", mva3);

        // iso deposits
        //Vetos 2011
        reco::isodeposit::AbsVetos vetos2010Charged;
        reco::isodeposit::AbsVetos vetos2010Neutral;
        reco::isodeposit::AbsVetos vetos2010Photons;
        reco::isodeposit::AbsVetos vetos2011Charged;
        reco::isodeposit::AbsVetos vetos2011Neutral;
        reco::isodeposit::AbsVetos vetos2011Photons;
        reco::isodeposit::AbsVetos vetos2011EECharged;
        reco::isodeposit::AbsVetos vetos2011EENeutral;
        reco::isodeposit::AbsVetos vetos2011EEPhotons;

        vetos2010Charged.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.01));
        vetos2010Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));
        vetos2010Neutral.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.08));
        vetos2010Neutral.push_back(new reco::isodeposit::ThresholdVeto(1.0));
        vetos2010Photons.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.05));
        vetos2010Photons.push_back(new reco::isodeposit::ThresholdVeto(1.0));

        vetos2011Charged.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.03)); //0.01
        vetos2011Charged.push_back(new reco::isodeposit::ThresholdVeto(0.5));//0.0
        vetos2011Neutral.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.08)); //0.01
        vetos2011Neutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));
        vetos2011Photons.push_back(new reco::isodeposit::ConeVeto(
                                       reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.05)); //0.01
        vetos2011Photons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

        vetos2011EECharged.push_back(new reco::isodeposit::ConeVeto(
                                         reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.03)); //0.015
        vetos2011EECharged.push_back(new reco::isodeposit::ThresholdVeto(0.5));//0.0
        vetos2011EENeutral.push_back(new reco::isodeposit::ConeVeto(
                                         reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.08)); //0.01
        vetos2011EENeutral.push_back(new reco::isodeposit::ThresholdVeto(0.5));
        vetos2011EEPhotons.push_back(new reco::isodeposit::ConeVeto(
                                         reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.05)); //0.08
        vetos2011EEPhotons.push_back(new reco::isodeposit::ThresholdVeto(0.5));

        //Vetos 2012
        reco::isodeposit::AbsVetos vetos2012EBPFIdCharged;
        reco::isodeposit::AbsVetos vetos2012EBPFIdNeutral;
        reco::isodeposit::AbsVetos vetos2012EBPFIdPhotons;
        reco::isodeposit::AbsVetos vetos2012EEPFIdCharged;
        reco::isodeposit::AbsVetos vetos2012EEPFIdNeutral;
        reco::isodeposit::AbsVetos vetos2012EEPFIdPhotons;

        reco::isodeposit::AbsVetos vetos2012EBNoPFIdCharged;
        reco::isodeposit::AbsVetos vetos2012EBNoPFIdNeutral;
        reco::isodeposit::AbsVetos vetos2012EBNoPFIdPhotons;
        reco::isodeposit::AbsVetos vetos2012EENoPFIdCharged;
        reco::isodeposit::AbsVetos vetos2012EENoPFIdNeutral;
        reco::isodeposit::AbsVetos vetos2012EENoPFIdPhotons;

        // safe recommended: PFId
        vetos2012EBPFIdCharged.push_back(new reco::isodeposit::ConeVeto(
                                             reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.010));
        vetos2012EBPFIdPhotons.push_back(new reco::isodeposit::ConeVeto(
                                             reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.08));
        vetos2012EEPFIdCharged.push_back(new reco::isodeposit::ConeVeto(
                                             reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.015));
        vetos2012EEPFIdPhotons.push_back(new reco::isodeposit::ConeVeto(
                                             reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.08));

        // POG recommended: NoPFId
        vetos2012EENoPFIdCharged.push_back(new reco::isodeposit::ConeVeto(
                                               reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.015));
        vetos2012EENoPFIdPhotons.push_back(new reco::isodeposit::ConeVeto(
                                               reco::isodeposit::Direction(aElectron.eta(), aElectron.phi()), 0.08));

        //2011
        float chIso03v1 =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2010Charged).first;
        float nhIso03v1 =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2010Neutral).first;
        float phIso03v1 =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2010Photons).first;
        float nhIsoPU03v1 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2010Neutral).first;
        float phIsoPU03v1 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2010Photons).first;

        float chIso04v1 =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2010Charged).first;
        float nhIso04v1 =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2010Neutral).first;
        float phIso04v1 =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2010Photons).first;
        float nhIsoPU04v1 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2010Neutral).first;
        float phIsoPU04v1 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2010Photons).first;

        float chIso03v2 =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2011Charged).first;
        float nhIso03v2 =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2011Neutral).first;
        float phIso03v2 =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2011Photons).first;
        float nhIsoPU03v2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2011Neutral).first;
        float phIsoPU03v2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2011Photons).first;

        float chIso04v2 =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2011Charged).first;
        float nhIso04v2 =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2011Neutral).first;
        float phIso04v2 =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2011Photons).first;
        float nhIsoPU04v2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011Neutral).first;
        float phIsoPU04v2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011Photons).first;

        float chIso03EEv2 =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2011EECharged).first;
        float nhIso03EEv2 =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2011EENeutral).first;
        float phIso03EEv2 =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2011EEPhotons).first;
        float nhIsoPU03EEv2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2011EENeutral).first;
        float phIsoPU03EEv2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2011EEPhotons).first;

        float chIso04EEv2 =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2011EECharged).first;
        float nhIso04EEv2 =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2011EENeutral).first;
        float phIso04EEv2 =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2011EEPhotons).first;
        float nhIsoPU04EEv2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011EENeutral).first;
        float phIsoPU04EEv2 =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2011EEPhotons).first;

        //2012
        float chIso03EBPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2012EBPFIdCharged).first;
        float nhIso03EBPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2012EBPFIdNeutral).first;
        float phIso03EBPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2012EBPFIdPhotons).first;
        float nhIsoPU03EBPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EBPFIdNeutral).first;
        //float phIsoPU03EBPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EBPFIdPhotons).first;

        float chIso03EEPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2012EEPFIdCharged).first;
        float nhIso03EEPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2012EEPFIdNeutral).first;
        float phIso03EEPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2012EEPFIdPhotons).first;
        float nhIsoPU03EEPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EEPFIdNeutral).first;
        //float phIsoPU03EEPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EEPFIdPhotons).first;

        float chIso03EBNoPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2012EBNoPFIdCharged).first;
        float nhIso03EBNoPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2012EBNoPFIdNeutral).first;
        float phIso03EBNoPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2012EBNoPFIdPhotons).first;
        float nhIsoPU03EBNoPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EBNoPFIdNeutral).first;
        //float phIsoPU03EBNoPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EBNoPFIdPhotons).first;

        float chIso03EENoPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.3, vetos2012EENoPFIdCharged).first;
        float nhIso03EENoPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.3, vetos2012EENoPFIdNeutral).first;
        float phIso03EENoPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.3, vetos2012EENoPFIdPhotons).first;
        float nhIsoPU03EENoPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EENoPFIdNeutral).first;
        //float phIsoPU03EENoPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.3, vetos2012EENoPFIdPhotons).first;

        //////////////////////////////////////////////////////////////

        float AllChIso04EBPFId =
            aElectron.isoDeposit(pat::User1Iso)->depositAndCountWithin(0.4, vetos2012EBPFIdCharged).first;
        float chIso04EBPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2012EBPFIdCharged).first;
        float nhIso04EBPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2012EBPFIdNeutral).first;
        float phIso04EBPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2012EBPFIdPhotons).first;
        float nhIsoPU04EBPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EBPFIdNeutral).first;
        //float phIsoPU04EBPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EBPFIdPhotons).first;

        float AllChIso04EEPFId =
            aElectron.isoDeposit(pat::User1Iso)->depositAndCountWithin(0.4, vetos2012EEPFIdCharged).first;
        float chIso04EEPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2012EEPFIdCharged).first;
        float nhIso04EEPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2012EEPFIdNeutral).first;
        float phIso04EEPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2012EEPFIdPhotons).first;
        float nhIsoPU04EEPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EEPFIdNeutral).first;
        //float phIsoPU04EEPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EEPFIdPhotons).first;

        float AllChIso04EBNoPFId =
            aElectron.isoDeposit(pat::User1Iso)->depositAndCountWithin(0.4, vetos2012EBNoPFIdCharged).first;
        float chIso04EBNoPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2012EBNoPFIdCharged).first;
        float nhIso04EBNoPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2012EBNoPFIdNeutral).first;
        float phIso04EBNoPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2012EBNoPFIdPhotons).first;
        float nhIsoPU04EBNoPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EBNoPFIdNeutral).first;
        //float phIsoPU04EBNoPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EBNoPFIdPhotons).first;

        float AllChIso04EENoPFId =
            aElectron.isoDeposit(pat::User1Iso)->depositAndCountWithin(0.4, vetos2012EENoPFIdCharged).first;
        float chIso04EENoPFId =
            aElectron.isoDeposit(pat::PfChargedHadronIso)->depositAndCountWithin(0.4, vetos2012EENoPFIdCharged).first;
        float nhIso04EENoPFId =
            aElectron.isoDeposit(pat::PfNeutralHadronIso)->depositAndCountWithin(0.4, vetos2012EENoPFIdNeutral).first;
        float phIso04EENoPFId =
            aElectron.isoDeposit(pat::PfGammaIso)->depositAndCountWithin(0.4, vetos2012EENoPFIdPhotons).first;
        float nhIsoPU04EENoPFId =
            aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EENoPFIdNeutral).first;
        //float phIsoPU04EENoPFId =
        //  aElectron.isoDeposit(pat::PfAllParticleIso)->depositAndCountWithin(0.4, vetos2012EENoPFIdPhotons).first;

        //////////////////////////////////////////////////////////////////////////////////////////////

        float chIso03 = pfId ?
                        (aElectron.isEB()) * chIso03EBPFId + (aElectron.isEE()) * chIso03EEPFId :
                        (aElectron.isEB()) * chIso03EBNoPFId + (aElectron.isEE()) * chIso03EENoPFId ;
        float nhIso03 = pfId ?
                        (aElectron.isEB()) * nhIso03EBPFId + (aElectron.isEE()) * nhIso03EEPFId :
                        (aElectron.isEB()) * nhIso03EBNoPFId + (aElectron.isEE()) * nhIso03EENoPFId ;
        float phIso03 = pfId ?
                        (aElectron.isEB()) * phIso03EBPFId + (aElectron.isEE()) * phIso03EEPFId :
                        (aElectron.isEB()) * phIso03EBNoPFId + (aElectron.isEE()) * phIso03EENoPFId ;
        float nhIsoPU03 = pfId ?
                          (aElectron.isEB()) * nhIsoPU03EBPFId + (aElectron.isEE()) * nhIsoPU03EEPFId :
                          (aElectron.isEB()) * nhIsoPU03EBNoPFId + (aElectron.isEE()) * nhIsoPU03EENoPFId ;

        float AllChIso04 = pfId ?
                        (aElectron.isEB()) * AllChIso04EBPFId + (aElectron.isEE()) * AllChIso04EEPFId :
                        (aElectron.isEB()) * AllChIso04EBNoPFId + (aElectron.isEE()) * AllChIso04EENoPFId ;
        float chIso04 = pfId ?
                        (aElectron.isEB()) * chIso04EBPFId + (aElectron.isEE()) * chIso04EEPFId :
                        (aElectron.isEB()) * chIso04EBNoPFId + (aElectron.isEE()) * chIso04EENoPFId ;
        float nhIso04 = pfId ?
                        (aElectron.isEB()) * nhIso04EBPFId + (aElectron.isEE()) * nhIso04EEPFId :
                        (aElectron.isEB()) * nhIso04EBNoPFId + (aElectron.isEE()) * nhIso04EENoPFId ;
        float phIso04 = pfId ?
                        (aElectron.isEB()) * phIso04EBPFId + (aElectron.isEE()) * phIso04EEPFId :
                        (aElectron.isEB()) * phIso04EBNoPFId + (aElectron.isEE()) * phIso04EENoPFId ;
        float nhIsoPU04 = pfId ?
                          (aElectron.isEB()) * nhIsoPU04EBPFId + (aElectron.isEE()) * nhIsoPU04EEPFId :
                          (aElectron.isEB()) * nhIsoPU04EBNoPFId + (aElectron.isEE()) * nhIsoPU04EENoPFId ;

        //2011
        aElectron.addUserFloat("PFRelIso04v1", (chIso04v1 + nhIso04v1 + phIso04v1) / aElectron.pt());
        aElectron.addUserFloat("PFRelIso03v1", (chIso03v1 + nhIso03v1 + phIso03v1) / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB04v1", (chIso04v1 + std::max(nhIso04v1 + phIso04v1 - 0.5 * 0.5 *
                                                                       (nhIsoPU04v1 + phIsoPU04v1), 0.0))
                                                 / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB03v1", (chIso03v1 + std::max(nhIso03v1 + phIso03v1 - 0.5 * 0.5 *
                                                                       (nhIsoPU03v1 + phIsoPU03v1), 0.0))
                                                 / aElectron.pt());

        aElectron.addUserFloat("PFRelIso04v2", (chIso04v2 + nhIso04v2 + phIso04v2) / aElectron.pt());
        aElectron.addUserFloat("PFRelIso03v2", (chIso03v2 + nhIso03v2 + phIso03v2) / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB04v2", (chIso04v2 + std::max(nhIso04v2 + phIso04v2 - 0.5 * 0.5 *
                                                                       (nhIsoPU04v2 + phIsoPU04v2), 0.0))
                                                 / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB03v2", (chIso03v2 + std::max(nhIso03v2 + phIso03v2 - 0.5 * 0.5 *
                                                                       (nhIsoPU03v2 + phIsoPU03v2), 0.0))
                                                 / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB04v3", (aElectron.isEB()) * (chIso04v2 +
                                                                       std::max(nhIso04v2 + phIso04v2 - 0.5 * 0.5 *
                                                                                (nhIsoPU04v2 + phIsoPU04v2), 0.0))
                                                 / aElectron.pt() +
                                                 (aElectron.isEE()) * (chIso04EEv2 +
                                                                       std::max(nhIso04EEv2 + phIso04EEv2 - 0.5 * 0.5 *
                                                                                (nhIsoPU04EEv2 + phIsoPU04EEv2), 0.0))
                                                 / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB03v3", (int(aElectron.isEB())) * (chIso03v2 +
                                                                            std::max(nhIso03v2 + phIso03v2 - 0.5 * 0.5 *
                                                                                     (nhIsoPU03v2 + phIsoPU03v2), 0.0))
                                                 / aElectron.pt() +
                                                 (int(aElectron.isEE())) * (chIso03EEv2 +
                                                                            std::max(nhIso03EEv2 + phIso03EEv2 - 0.5 *
                                                                                     0.5 * (nhIsoPU03EEv2 +
                                                                                            phIsoPU03EEv2), 0.0))
                                                 / aElectron.pt());

        //2012
        aElectron.addUserFloat("PFRelIso03",
                               (chIso03 + nhIso03 + phIso03) / aElectron.pt());
        aElectron.addUserFloat("PFRelIsoDB03",
                               (chIso03 + std::max(nhIso03 + phIso03 - 0.5 * (nhIsoPU03), 0.0)) / aElectron.pt());

        aElectron.addUserFloat("PFRelIso04",
                               (chIso04 + nhIso04 + phIso04) / aElectron.pt());
        // with charged hadrons only
//        aElectron.addUserFloat("PFRelIsoDB04",
//                               (chIso04 + std::max(nhIso04 + phIso04 - 0.5 * (nhIsoPU04), 0.0)) / aElectron.pt());
        // legacy definition
        aElectron.addUserFloat("PFRelIsoDB04",
                               (AllChIso04 + std::max(nhIso04 + phIso04 - 0.5 * (nhIsoPU04), 0.0)) / aElectron.pt());

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
        for(unsigned int i = 0; i < vetos2011EECharged.size(); i++) {
            delete vetos2011EECharged[i];
        }
        for(unsigned int i = 0; i < vetos2011EENeutral.size(); i++) {
            delete vetos2011EENeutral[i];
            delete vetos2011EEPhotons[i];
        }

        for(unsigned int i = 0; i < vetos2012EBPFIdCharged.size(); i++) {
            delete vetos2012EBPFIdCharged[i];
        }
        for(unsigned int i = 0; i < vetos2012EBPFIdNeutral.size(); i++) {
            delete vetos2012EBPFIdNeutral[i];
        }
        for(unsigned int i = 0; i < vetos2012EBPFIdPhotons.size(); i++) {
            delete vetos2012EBPFIdPhotons[i];
        }
        for(unsigned int i = 0; i < vetos2012EEPFIdCharged.size(); i++) {
            delete vetos2012EEPFIdCharged[i];
        }
        for(unsigned int i = 0; i < vetos2012EEPFIdNeutral.size(); i++) {
            delete vetos2012EEPFIdNeutral[i];
        }
        for(unsigned int i = 0; i < vetos2012EEPFIdPhotons.size(); i++) {
            delete vetos2012EEPFIdPhotons[i];
        }
        for(unsigned int i = 0; i < vetos2012EBNoPFIdCharged.size(); i++) {
            delete vetos2012EBNoPFIdCharged[i];
        }
        for(unsigned int i = 0; i < vetos2012EBNoPFIdNeutral.size(); i++) {
            delete vetos2012EBNoPFIdNeutral[i];
        }
        for(unsigned int i = 0; i < vetos2012EBNoPFIdPhotons.size(); i++) {
            delete vetos2012EBNoPFIdPhotons[i];
        }
        for(unsigned int i = 0; i < vetos2012EENoPFIdCharged.size(); i++) {
            delete vetos2012EENoPFIdCharged[i];
        }
        for(unsigned int i = 0; i < vetos2012EENoPFIdNeutral.size(); i++) {
            delete vetos2012EENoPFIdNeutral[i];
        }
        for(unsigned int i = 0; i < vetos2012EENoPFIdPhotons.size(); i++) {
            delete vetos2012EENoPFIdPhotons[i];
        }

        aElectron.addUserInt("pfId", pfId);
        aElectron.addUserFloat("isInRun", iEvent.run());

        electronsUserEmbeddedColl->push_back(aElectron);

    }

    iEvent.put( electronsUserEmbeddedColl );
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronsUserEmbedder);
