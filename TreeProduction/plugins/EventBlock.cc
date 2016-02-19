/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


#include "h-tautau/PatProduction/interface/PatVertex.h"
#include "h-tautau/TreeProduction/interface/Event.h"

class EventBlock : public edm::EDAnalyzer {
public:
    explicit EventBlock(const edm::ParameterSet& iConfig) :
        _l1InputTag(iConfig.getParameter<edm::InputTag>("l1InputTag")),
        _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexInputTag")),
        _vtxMinNDOF(iConfig.getParameter<unsigned int>("vertexMinimumNDOF")),
        _vtxMaxAbsZ(iConfig.getParameter<double>("vertexMaxAbsZ")),
        _vtxMaxd0(iConfig.getParameter<double>("vertexMaxd0")),
        _hpTrackThreshold(iConfig.getParameter<double>("hpTrackThreshold")),
        eventTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { eventTree.Write(); }
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);

private:
    const edm::InputTag   _l1InputTag;
    const edm::InputTag   _vtxInputTag;
    const unsigned int    _vtxMinNDOF;
    const double          _vtxMaxAbsZ, _vtxMaxd0;
    const double          _hpTrackThreshold;

    ntuple::EventTree eventTree;
};

void EventBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)
{
    eventTree.run()   = iEvent.id().run();
    eventTree.EventId() = iEvent.id().event();
    eventTree.lumis() = iEvent.id().luminosityBlock();
    eventTree.bunch() = iEvent.bunchCrossing();
    eventTree.orbit() = iEvent.orbitNumber();

    eventTree.unixTime() = iEvent.time().unixTime();
    eventTree.microsecondOffset() = iEvent.time().microsecondOffset();
    eventTree.isdata() = iEvent.isRealData();

    edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
    iEvent.getByLabel(_l1InputTag, l1GtReadoutRecord);

    // Technical Trigger Part
    if (l1GtReadoutRecord.isValid()) {
        edm::LogInfo("EventBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: "
                                   << _l1InputTag;

        L1GtFdlWord fdlWord = l1GtReadoutRecord->gtFdlWord();
        eventTree.isPhysDeclared() = fdlWord.physicsDeclared() == 1;


        // BPTX0
        eventTree.isBPTX0() = l1GtReadoutRecord->technicalTriggerWord()[0];


        // MinBias
        eventTree.isBSCMinBias() = l1GtReadoutRecord->technicalTriggerWord()[40] ||
                                   l1GtReadoutRecord->technicalTriggerWord()[41];


        // BeamHalo
        eventTree.isBSCBeamHalo() =
            (l1GtReadoutRecord->technicalTriggerWord()[36] || l1GtReadoutRecord->technicalTriggerWord()[37] ||
             l1GtReadoutRecord->technicalTriggerWord()[38] || l1GtReadoutRecord->technicalTriggerWord()[39]) ||
            ((l1GtReadoutRecord->technicalTriggerWord()[42] && !l1GtReadoutRecord->technicalTriggerWord()[43]) ||
             (l1GtReadoutRecord->technicalTriggerWord()[43] && !l1GtReadoutRecord->technicalTriggerWord()[42])) ;

    } else {
        edm::LogError("EventBlock") << "Error >> Failed to get L1GlobalTriggerReadoutRecord for label:"
                                    << _l1InputTag;
        throw std::runtime_error("Failed to get L1GlobalTriggerReadoutRecord for label.");
    }

    // Good Primary Vertex Part
    edm::Handle<pat::VertexCollection> primaryVertices;
    iEvent.getByLabel(_vtxInputTag, primaryVertices);

    eventTree.isPrimaryVertex() = false;
    if (primaryVertices.isValid()) {
        edm::LogInfo("EventBlock") << "Total # Primary Vertices: " << primaryVertices->size();
        for (auto it = primaryVertices->begin(); it != primaryVertices->end(); ++it) {
            if (!(it->isFake()) && it->ndof() > _vtxMinNDOF &&
                    fabs(it->z()) <= _vtxMaxAbsZ && fabs(it->position().rho()) <= _vtxMaxd0
               ) {
                eventTree.isPrimaryVertex() = true;
                break;
            }
        }
    } else {
        edm::LogError("EventBlock") << "Error >> Failed to get VertexCollection for label:"
                                    << _vtxInputTag;
        throw std::runtime_error("Failed to get VertexCollection for label.");
    }


    // Access PU information
    if (!iEvent.isRealData()) {
        edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
        iEvent.getByLabel(edm::InputTag("addPileupInfo","","HLT"), PupInfo);

        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            eventTree.bunchCrossing().push_back(PVI->getBunchCrossing());
            eventTree.nPU().push_back(PVI->getPU_NumInteractions());
            eventTree.trueNInt().push_back(PVI->getTrueNumInteractions());
        }

        // More info about PU is here:
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
    }
    eventTree.Fill();

}

DEFINE_FWK_MODULE(EventBlock);
