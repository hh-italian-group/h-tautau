/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "h-tautau/TreeProduction/interface/MET.h"

class METBlock : public edm::EDAnalyzer {
public:
    typedef std::shared_ptr<ntuple::METTree> TreePtr;
    typedef std::pair<edm::InputTag, TreePtr> TreeDescriptor;
    typedef std::vector<TreeDescriptor> TreeDescriptorVector;

    explicit METBlock(const edm::ParameterSet& iConfig)
    {
        const edm::ParameterSet& sources = iConfig.getParameterSet("metSrc");
        const auto names = sources.getParameterNames();
        edm::Service<TFileService> fileService;
        for(const std::string& name : names) {
            const edm::InputTag tag = sources.getParameter<edm::InputTag>(name);
            const TreeDescriptor descriptor(tag, TreePtr(new ntuple::METTree(name, &fileService->file(), false)));
            descriptors.push_back(descriptor);
        }
    }

private:
    virtual void endJob()
    {
        for(const auto& descriptor : descriptors)
            descriptor.second->Write();
    }

    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    TreeDescriptorVector descriptors;
};

void METBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    for(const auto& descriptor : descriptors) {
        const edm::InputTag& inputTag = descriptor.first;
        ntuple::METTree& METTree = *descriptor.second.get();
        METTree.RunId() = iEvent.id().run();
        METTree.LumiBlock() = iEvent.id().luminosityBlock();
        METTree.EventId() = iEvent.id().event();

        edm::Handle<std::vector<pat::MET> > mets;
        iEvent.getByLabel(inputTag, mets);

        if (!mets.isValid()) {
            edm::LogError("METBlock") << "Error: Failed to get METCollection for label: " << inputTag;
            throw std::runtime_error("Failed to get METCollection.");
        }
        edm::LogInfo("METBlock") << "Total # METs: " << mets->size();
        for (const pat::MET& MET : *mets) {
            METTree.pt() = MET.pt();
            METTree.phi() = MET.phi();
            METTree.sumEt()  = MET.sumEt();
           //            METTree.pt_uncorrected() = MET.uncorrectedPt(pat::MET::uncorrALL);
            //            METTree.phi_uncorrected() = MET.uncorrectedPhi(pat::MET::uncorrALL);
            //            METTree.sumEt_uncorrected() = MET.sumEt() - MET.corSumEt(pat::MET::uncorrALL);
               //Compatibility with 74X
            METTree.pt_uncorrected() = MET.uncorPt();
            METTree.phi_uncorrected() = MET.uncorPhi();
            METTree.sumEt_uncorrected() = MET.sumEt() - MET.uncorSumEt();
            METTree.significanceMatrix() = ntuple::SignificanceMatrixToVector(MET.getSignificanceMatrix());
            METTree.Fill();
            METTree.Fill();
        }
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
