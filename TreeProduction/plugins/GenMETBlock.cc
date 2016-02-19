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
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"

#include "h-tautau/TreeProduction/interface/GenMET.h"

class GenMETBlock : public edm::EDAnalyzer {
public:
    explicit GenMETBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("genMETSrc")),
        genMETTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { genMETTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
  int _verbosity;
  edm::InputTag _inputTag;
  ntuple::GenMETTree genMETTree;
};

void GenMETBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    genMETTree.RunId() = iEvent.id().run();
    genMETTree.LumiBlock() = iEvent.id().luminosityBlock();
    genMETTree.EventId() = iEvent.id().event();

    edm::Handle<reco::GenMETCollection> mets;
    iEvent.getByLabel(_inputTag, mets);
    if (!mets.isValid()) {
        edm::LogError("GenMETBlock") << "Error >>  Failed to get GenMETCollection for label: "
                                     << _inputTag;
        throw std::runtime_error("Failed to get GenMETCollection.");
    }
      edm::LogInfo("GenMETBlock") << "Total # GenMETs: " << mets->size();
      for (const reco::GenMET& genMET : *mets) {
        genMETTree.met()    = genMET.pt();
        genMETTree.metphi() = genMET.phi();
        genMETTree.sumet()  = genMET.sumEt();

        genMETTree.Fill();
    } 
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenMETBlock);
