/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TPRegexp.h"

#include "h-tautau/TreeProduction/interface/Trigger.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include <string>
#include <vector>
#include <limits>

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

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

static const unsigned int NmaxL1AlgoBit = 128;
static const unsigned int NmaxL1TechBit = 64;

class TriggerBlock : public edm::EDAnalyzer
{
public:
    // Constructor
    TriggerBlock(const edm::ParameterSet& iConfig) :
      _verbosity(iConfig.getParameter<int>("verbosity")),
      _l1InputTag(iConfig.getParameter<edm::InputTag>("l1InputTag")),
      _hltInputTag(iConfig.getParameter<edm::InputTag>("hltInputTag")),
      _hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("hltPathsOfInterest")),
      triggerTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void endJob() { triggerTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);


private:
  int _verbosity;
  const edm::InputTag _l1InputTag;
  const edm::InputTag _hltInputTag;
  const std::vector<std::string> _hltPathsOfInterest;
  HLTConfigProvider hltConfig;
  ntuple::TriggerTree triggerTree;


};


void TriggerBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig.init(iRun, iSetup, _hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerBlock") << "HLT config with process name "
                                 << _hltInputTag.process() << " successfully extracted";
  }
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerBlock") << "Error! HLT config extraction with process name "
                                  << _hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
    throw std::runtime_error("Failed to get Trigger Collection");
  }

}


void TriggerBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    triggerTree.RunId() = iEvent.id().run();
    triggerTree.LumiBlock() = iEvent.id().luminosityBlock();
    triggerTree.EventId() = iEvent.id().event();


  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(_l1InputTag, l1GtReadoutRecord);

    if (!l1GtReadoutRecord.isValid()) {
      edm::LogError("TriggerBlock") << "Error >> Failed to get L1GlobalTriggerReadoutRecord for label: "
                                    << _l1InputTag;
      throw std::runtime_error("Failed to get L1Trigger Collection");
    }
    edm::LogInfo("TriggerBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: "
                                 << _l1InputTag;
    triggerTree.l1physbits() = l1GtReadoutRecord->decisionWord();
    triggerTree.l1techbits() = l1GtReadoutRecord->technicalTriggerWord();

    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(_hltInputTag, triggerResults);

    if (!triggerResults.isValid()) {
        edm::LogError("TriggerBlock") << "Error >> Failed to get TriggerResults for label: "
                                    << _hltInputTag;
        throw std::runtime_error("Failed to get TriggerResult Collection");
    }

    edm::LogInfo("TriggerBlock") << "Successfully obtained " << _hltInputTag;
    const std::vector<std::string>& pathList = hltConfig.triggerNames();
    for (std::vector<std::string>::const_iterator it = pathList.begin();
                                                 it != pathList.end(); ++it) {
        if (_hltPathsOfInterest.size()) {
            int nmatch = 0;
            for (std::vector<std::string>::const_iterator kt = _hltPathsOfInterest.begin();
                                                     kt != _hltPathsOfInterest.end(); ++kt) {
                nmatch += TPRegexp(*kt).Match(*it);
        }
            if (!nmatch) continue;
        }
        triggerTree.hltpaths().push_back(*it);

        bool fired = false;
        unsigned int index = hltConfig.triggerIndex(*it);
        if (index < triggerResults->size()) {
            if (triggerResults->accept(index)) fired = true;
        }
        else {
            edm::LogInfo("TriggerBlock") << "Requested HLT path \"" << (*it) << "\" does not exist";
        }
        triggerTree.hltresults().push_back(fired);

        unsigned prescale = std::numeric_limits<unsigned>::max();
        if (hltConfig.prescaleSet(iEvent, iSetup) < 0) {
            edm::LogError("TriggerBlock") << "Error >> The prescale set index number could not be obtained";
        }
        else {
            prescale = hltConfig.prescaleValue(iEvent, iSetup, *it);
        }
        triggerTree.hltprescales().push_back(prescale);

        if (_verbosity)
        std::cout << ">>> Path: " << (*it)
                  << ", prescale: " << prescale
                  << ", fired: " << fired
                  << std::endl;
    }
    triggerTree.Fill();

}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerBlock);
