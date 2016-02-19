/*!
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include <iostream>
#include <algorithm>

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

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
 
#include "h-tautau/TreeProduction/interface/Jet.h"
#include "h-tautau/TreeProduction/interface/TriggerTools.h"

#define SIMPLE_VAR(type, name) jetTree.name() = jet.bDiscriminator(#name);

namespace {
PFJetIDSelectionFunctor pfjetIDLoose(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE);
PFJetIDSelectionFunctor pfjetIDTight(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
pat::strbitset retpf = pfjetIDLoose.getBitTemplate();
}

class JetBlock : public edm::EDAnalyzer {
public:
    explicit JetBlock(const edm::ParameterSet& iConfig) :
        _inputTag(iConfig.getParameter<edm::InputTag>("jetSrc")),
        jetTree(&edm::Service<TFileService>()->file(), false) {}

private:
    virtual void endJob() { jetTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    edm::InputTag _inputTag;
    ntuple::JetTree jetTree;
};

void JetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    jetTree.RunId() = iEvent.id().run();
    jetTree.LumiBlock() = iEvent.id().luminosityBlock();
    jetTree.EventId() = iEvent.id().event();

    edm::Handle<edm::View<pat::Jet> > jets;
    iEvent.getByLabel(_inputTag, jets);
    if (!jets.isValid()) {
        edm::LogError("JetBlock") << "Error >> Failed to get pat::Jet collection for label: "
                                  << _inputTag;
        throw std::runtime_error("Failed to get pat::Jet collection");
    }

    unsigned int njets = jets->size();
    edm::LogInfo("JetBlock") << "Total # PAT Jets: " << njets;
    for (size_t i = 0; i < njets; ++i) {
      const pat::Jet& jet = jets->at(i);

      retpf.set(false);
      jetTree.passLooseID() = pfjetIDLoose(jet, retpf);

      retpf.set(false);
      jetTree.passTightID() = pfjetIDTight(jet, retpf);

      // fill in all the vectors
      jetTree.eta()        = jet.eta();
      jetTree.phi()        = jet.phi();
      jetTree.pt()         = jet.pt();
      jetTree.mass()     = jet.mass();
      jetTree.pt_raw()     = jet.correctedJet("Uncorrected").pt();
      jetTree.energy_raw() = jet.correctedJet("Uncorrected").energy();
      jetTree.partonFlavour() = jet.partonFlavour();

      // Jet identification in high pile-up environment
      jetTree.puIdMVA() = jet.userFloat("pileupJetIdProducer:fullDiscriminant");
      jetTree.puIdBits() = jet.userInt("pileupJetIdProducer:fullId"); // Bits: 0:Tight,1:Medium,2:Loose
      jetTree.puIdMVA_met() = jet.userFloat("pileupJetIdProducer:metDiscriminant");
      jetTree.puIdBits_met() = jet.userInt("pileupJetIdProducer:metId"); // Bits: 0:Tight,1:Medium,2:Loose
      jetTree.correction() = jet.userFloat("correction");

      jetTree.chargedEmEnergyFraction()     = jet.chargedEmEnergyFraction();
      jetTree.chargedHadronEnergyFraction() = jet.chargedHadronEnergyFraction();
      jetTree.chargedMuEnergyFraction()     = jet.chargedMuEnergyFraction();
      jetTree.electronEnergyFraction()      = jet.electronEnergy()/jet.energy();
      jetTree.muonEnergyFraction()          = jet.muonEnergyFraction();
      jetTree.neutralEmEnergyFraction()     = jet.neutralEmEnergyFraction();
      jetTree.neutralHadronEnergyFraction() = jet.neutralHadronEnergyFraction();
      jetTree.photonEnergyFraction()        = jet.photonEnergyFraction();

      jetTree.chargedHadronMultiplicity()   = jet.chargedHadronMultiplicity();
      jetTree.chargedMultiplicity()         = jet.chargedMultiplicity();
      jetTree.electronMultiplicity()        = jet.electronMultiplicity();
      jetTree.muonMultiplicity()            = jet.muonMultiplicity();
      jetTree.neutralHadronMultiplicity()   = jet.neutralHadronMultiplicity();
      jetTree.neutralMultiplicity()         = jet.neutralMultiplicity();
      jetTree.photonMultiplicity()          = jet.photonMultiplicity();

      jetTree.nConstituents()               = jet.numberOfDaughters();

      jetTree.matchedTriggerPaths() = CollectMatchedTriggerPaths(jet);

      B_TAG_DATA()

      jetTree.Fill();
    }
}

#undef SIMPLE_VAR
#undef B_TAG_DATA

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(JetBlock);
