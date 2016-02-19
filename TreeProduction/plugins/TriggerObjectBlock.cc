/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>
#include <algorithm>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerFilter.h"

#include "TMath.h"
#include "TPRegexp.h"

#include "h-tautau/TreeProduction/interface/TriggerObject.h"

class TriggerObjectBlock : public edm::EDAnalyzer {
public:
    explicit TriggerObjectBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _hltInputTag(iConfig.getParameter<edm::InputTag>("hltInputTag")),
        _triggerEventTag(iConfig.getParameter<edm::InputTag>("triggerEventTag")),
        _hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("hltPathsOfInterest")),
        _may10ReRecoData(iConfig.getParameter<bool>("May10ReRecoData")),
        _firingFlag(_may10ReRecoData),
        triggerObjectTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { triggerObjectTree.Write(); }
    virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    const edm::InputTag _hltInputTag;
    edm::InputTag _triggerEventTag;
    const std::string  _tagPathLabel;
    const std::string  _probePathLabel;
    const std::vector<std::string> _hltPathsOfInterest;
    bool _may10ReRecoData;
    bool _firingFlag;

    ntuple::TriggerObjectTree triggerObjectTree;
    HLTConfigProvider hltConfig;
};

void TriggerObjectBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    bool changed = true;
    if (hltConfig.init(iRun, iSetup, _hltInputTag.process(), changed)) {
        // if init returns TRUE, initialisation has succeeded!
        edm::LogInfo("TriggerObjectBlock") << "HLT config with process name "
                                   << _hltInputTag.process()
                                   << " successfully extracted";
    }
    else {
        // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
        // with the file and/or code and needs to be investigated!
        edm::LogError("TriggerObjectBlock") << "Error! HLT config extraction with process name "
                                    << _hltInputTag.process() << " failed";
        // In this case, all access methods will return empty values!
        throw std::runtime_error("HLT config extraction failed.");
    }
}

void TriggerObjectBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    triggerObjectTree.RunId() = iEvent.id().run();
    triggerObjectTree.LumiBlock() = iEvent.id().luminosityBlock();
    triggerObjectTree.EventId() = iEvent.id().event();

    if (_verbosity) {
        std::cout << setiosflags(std::ios::fixed);
        std::cout << "Indx     Eta     Phi      Pt  Energy            =Trigger path list=" << std::endl;
    }

    // trigger event
    edm::Handle<pat::TriggerEvent> triggerEvent;
    iEvent.getByLabel(_triggerEventTag, triggerEvent);

    // get the trigger objects corresponding to the used matching (HLT muons) and
    // loop over selected trigger objects
    pat::TriggerObjectRefVector myObjects(triggerEvent->objectRefs());
    for (pat::TriggerObjectRefVector::const_iterator it  = myObjects.begin(); it != myObjects.end(); ++it) {
        pat::TriggerPathRefVector myPaths = triggerEvent->objectPaths((*it));
        std::map <std::string, bool> pathInfoMap;

        for (pat::TriggerPathRefVector::const_iterator ipath  = myPaths.begin(); ipath != myPaths.end(); ++ipath) {
            std::string name = (**ipath).name();
            for (const std::string& path_int : _hltPathsOfInterest) {
                if (name.find(path_int) == std::string::npos) continue;
                bool matched = true;
                // Get the filters and access the L3 filter (needed for May10ReReco data)
                if (_may10ReRecoData) {
                    matched = false;
                    pat::TriggerFilterRefVector filters( triggerEvent->pathFilters( name, _firingFlag) );
                    if ( filters.empty() ) continue;
                    pat::TriggerFilterRef lastFilter( filters.at( filters.size() - 1 ) );
                    if ( triggerEvent->objectInFilter( (*it), lastFilter->label() ) ) matched = true;
                }
                if (matched) {
                    bool val = false;
                    if (triggerEvent->path(name)->wasRun() && triggerEvent->path(name)->wasAccept()) val = true;
                    pathInfoMap.insert(std::pair<std::string, bool> (name, val));
                }
            }
        }

        if (!pathInfoMap.size()) continue;

        triggerObjectTree.eta()    = (**it).eta();
        triggerObjectTree.phi()    = (**it).phi();
        triggerObjectTree.pt()     = (**it).pt();
        triggerObjectTree.mass() = (**it).mass();
        triggerObjectTree.pdgId() = (**it).pdgId();

        for (const auto& imap : pathInfoMap) {
            triggerObjectTree.pathNames().push_back(imap.first);
            triggerObjectTree.pathValues().push_back(imap.second);
        }

        triggerObjectTree.Fill();
    }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectBlock);
