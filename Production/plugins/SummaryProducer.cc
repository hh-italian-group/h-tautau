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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"

class SummaryProducer : public edm::EDAnalyzer {
public:
    using clock = std::chrono::system_clock;
    using GenId = ntuple::GenId;
    using GenEventCountMap = ntuple::GenEventCountMap;

    SummaryProducer(const edm::ParameterSet& cfg) :
        start(clock::now()),
        isMC(cfg.getParameter<bool>("isMC")),
        lheEventProduct_token(mayConsume<LHEEventProduct>(cfg.getParameter<edm::InputTag>("lheEventProduct"))),
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
        static constexpr int b_quark = 5;
        static const std::set<int> quarks_and_gluons = { 1, 2, 3, 4, 5, 6, 21 };

        summaryTuple().numberOfProcessedEvents++;

        if(isMC) {
            edm::Handle<GenEventInfoProduct> genEvent;
            event.getByToken(genEvent_token, genEvent);
            summaryTuple().totalWeight += genEvent->weight();

            edm::Handle<LHEEventProduct> lheEventProduct;
            event.getByToken(lheEventProduct_token, lheEventProduct);
            if(lheEventProduct.isValid()) {
                const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
                const std::vector<lhef::HEPEUP::FiveVector>& lheParticles = lheEvent.PUP;
                size_t n_partons = 0, n_b_partons = 0;
                double HT = 0;
                for(size_t n = 0; n < lheParticles.size(); ++n) {
                    const int absPdgId = std::abs(lheEvent.IDUP[n]);
                    const int status = lheEvent.ISTUP[n];
                    if(status != 1 || !quarks_and_gluons.count(absPdgId)) continue;
                    ++n_partons;
                    if(absPdgId == b_quark) ++n_b_partons;
                    HT += std::sqrt(std::pow(lheParticles[n][0], 2) + std::pow(lheParticles[n][1], 2));
                }
                const size_t ht10_bin = HT / 10;
                const GenId genId(n_partons, n_b_partons, ht10_bin);
                ++genEventCountMap[genId];
            }
        } else {
            summaryTuple().totalWeight += 1;
        }
        if(!tauId_names.size()) {
            edm::Handle<std::vector<pat::Tau>> taus;
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
        for(const auto& count_entry : genEventCountMap) {
            summaryTuple().lhe_n_partons.push_back(count_entry.first.n_partons);
            summaryTuple().lhe_n_b_partons.push_back(count_entry.first.n_b_partons);
            summaryTuple().lhe_ht10_bin.push_back(count_entry.first.ht10_bin);
            summaryTuple().lhe_n_events.push_back(count_entry.second);
        }
        const auto stop = clock::now();
        summaryTuple().exeTime = std::chrono::duration_cast<std::chrono::seconds>(stop - start).count();
        summaryTuple.Fill();
        summaryTuple.Write();
    }

private:
    const clock::time_point start;
    const bool isMC;

    edm::EDGetTokenT<LHEEventProduct> lheEventProduct_token;
    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<std::vector<pat::Tau>> taus_token;

    ntuple::SummaryTuple summaryTuple;
    std::unordered_set<std::string> tauId_names;
    GenEventCountMap genEventCountMap;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SummaryProducer);
