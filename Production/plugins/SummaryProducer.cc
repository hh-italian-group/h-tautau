/*! Creates tuple with the production summary.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"

class SummaryProducer : public edm::EDAnalyzer {
public:
    SummaryProducer(const edm::ParameterSet& cfg) :
        isMC(cfg.getParameter<bool>("isMC")),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        taus_token(mayConsume<std::vector<pat::Tau>>(cfg.getParameter<edm::InputTag>("taus"))),
        summaryTuple("summary", &edm::Service<TFileService>()->file(), false)
    {
        summaryTuple().numberOfProcessedEvents = 0;
        summaryTuple().totalWeight = 0.;
    }

private:
    virtual void analyze(const edm::Event& event, const edm::EventSetup&) override
    {
        summaryTuple().numberOfProcessedEvents++;

        if(isMC) {
            event.getByToken(genEvent_token, genEvent);
            summaryTuple().totalWeight += genEvent->weight();
        } else {
            summaryTuple().totalWeight += 1;
        }
        if(!tauId_names.size()) {
            event.getByToken(taus_token, taus);
            for(const auto& tau : *taus) {
                for(const auto& id : tau.tauIDs()) {
                    tauId_names.insert(id.first);
                }
            }
        }
    }

    virtual void endJob() override
    {
        for(const std::string& name : tauId_names) {
            const auto key = analysis::tools::hash(name);
            summaryTuple().tauId_names.push_back(name);
            summaryTuple().tauId_keys.push_back(key);
        }
        summaryTuple.Fill();
        summaryTuple.Write();
    }

private:
    const bool isMC;
    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<pat::Tau>> taus_token;

    edm::Handle<GenEventInfoProduct> genEvent;
    edm::Handle<std::vector<pat::Tau>> taus;

    ntuple::SummaryTuple summaryTuple;
    std::unordered_set<std::string> tauId_names;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SummaryProducer);
