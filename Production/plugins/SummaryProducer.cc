/*! Creates tuple with the production summary.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TauIdResults.h"
#include "h-tautau/Production/interface/GenTruthTools.h"


class SummaryProducer : public edm::EDAnalyzer {
public:
    using clock = std::chrono::system_clock;
    using GenId = ntuple::GenId;
    using GenEventCountMap = ntuple::GenEventCountMap;
    using GenEventTypeCountMap = ntuple::GenEventTypeCountMap;
    using Channel = analysis::Channel;

    SummaryProducer(const edm::ParameterSet& cfg) :
        start(clock::now()),
        isMC(cfg.getParameter<bool>("isMC")),
        saveGenTopInfo(cfg.getParameter<bool>("saveGenTopInfo")),
        saveGenBosonInfo(cfg.getParameter<bool>("saveGenBosonInfo")),
        lheEventProduct_token(mayConsume<LHEEventProduct>(cfg.getParameter<edm::InputTag>("lheEventProduct"))),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        topGenEvent_token(mayConsume<TtGenEvent>(cfg.getParameter<edm::InputTag>("topGenEvent"))),
        genParticles_token(mayConsume<std::vector<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        summaryTuple(ntuple::CreateSummaryTuple("summary", &edm::Service<TFileService>()->file(), false,
                                                ntuple::TreeState::Full))
    {
        (*summaryTuple)().numberOfProcessedEvents = 0;
        if(isMC)
            expressTuple = std::shared_ptr<ntuple::ExpressTuple>(
                    new ntuple::ExpressTuple("all_events", &edm::Service<TFileService>()->file(), false));

    }

private:
    virtual void analyze(const edm::Event& event, const edm::EventSetup&) override
    {
        (*summaryTuple)().numberOfProcessedEvents++;

        if(!isMC)
            return;

        edm::Handle<GenEventInfoProduct> genEvent;
        event.getByToken(genEvent_token, genEvent);

        edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
        event.getByToken(puInfo_token, puInfo);
        (*expressTuple)().npu = analysis::gen_truth::GetNumberOfPileUpInteractions(puInfo);
        (*expressTuple)().run = event.id().run();
        (*expressTuple)().lumi = event.id().luminosityBlock();
        (*expressTuple)().evt = event.id().event();
        (*expressTuple)().genEventWeight = genEvent->weight();
        (*expressTuple)().genEventWeight = isMC ? genEvent->weight() : 1;
        (*expressTuple)().genEventLHEWeight = isMC && genEvent->weights().size() > 1 ? genEvent->weights()[1] : 1;

        if(genEvent->weights().size() == 14 || genEvent->weights().size() == 46) {
            const double nominal = genEvent->weights()[1];  // Called 'Baseline' in GenLumiInfoHeader
            for (size_t i = 6; i < 10; ++i)
                (*expressTuple)().genEventPSWeights.push_back(genEvent->weights()[i] / nominal);
        }

        (*expressTuple)().gen_top_pt = ntuple::DefaultFillValue<Float_t>();
        (*expressTuple)().gen_topBar_pt = ntuple::DefaultFillValue<Float_t>();
        (*expressTuple)().lhe_H_m = ntuple::DefaultFillValue<Float_t>();

        if(saveGenTopInfo || saveGenBosonInfo)
            FillGenParticleInfo(event);

        edm::Handle<LHEEventProduct> lheEventProduct;
        event.getByToken(lheEventProduct_token, lheEventProduct);
        if(lheEventProduct.isValid()) {
            const auto lheSummary = analysis::gen_truth::ExtractLheSummary(*lheEventProduct);
            const size_t ht10_bin = lheSummary.HT / 10;
            const GenId genId(lheSummary.n_partons, lheSummary.n_b_partons, ht10_bin);
            ++genEventCountMap[genId];
            (*expressTuple)().lhe_H_m = lheSummary.m_H;
            (*expressTuple)().lhe_hh_m = lheSummary.m_hh;
            (*expressTuple)().lhe_hh_cosTheta = lheSummary.cosTheta_hh;
            (*expressTuple)().lhe_n_partons = lheSummary.n_partons;
            (*expressTuple)().lhe_n_b_partons = lheSummary.n_b_partons;
            (*expressTuple)().lhe_ht10_bin = ht10_bin;
        }

        expressTuple->Fill();
    }

    virtual void endJob() override
    {
        const auto stop = clock::now();
        (*summaryTuple)().exeTime = std::chrono::duration_cast<std::chrono::seconds>(stop - start).count();
        summaryTuple->Fill();
        summaryTuple->Write();
        if(expressTuple)
            expressTuple->Write();
    }

    void FillGenParticleInfo(const edm::Event& event)
    {
        using analysis::GenEventType;
        static constexpr int electronPdgId = 11, muonPdgId = 13, tauPdgId = 15;
        static constexpr int electronNeutrinoPdgId = 12, muonNeutrinoPdgId = 14, tauNeutrinoPdgId = 16;
        static constexpr int topPdgId = 6;
        static const std::set<int> chargedLeptons = { electronPdgId, muonPdgId, tauPdgId };
        static const std::set<int> neutralLeptons = { electronNeutrinoPdgId, muonNeutrinoPdgId, tauNeutrinoPdgId };
        static const std::set<int> bosons = { 23, 24, 25, 35 };

        std::vector<const reco::GenParticle*> particles_to_store, mothers_to_store;

        const auto findStoredMother = [&](const reco::GenParticle& p) -> const reco::GenParticle* {
            for(auto iter = particles_to_store.rbegin(); iter != particles_to_store.rend(); ++iter) {
                if(analysis::gen_truth::CheckAncestry(p, **iter))
                    return *iter;
            }
            return nullptr;
        };

        edm::Handle<std::vector<reco::GenParticle>> genParticles;
        event.getByToken(genParticles_token, genParticles);
        if(genParticles.isValid()) {
            for(const auto& particle : *genParticles) {
                const auto& flag = particle.statusFlags();
                if(!flag.isPrompt() || !flag.isLastCopy() || !flag.fromHardProcess()) continue;
                const int abs_pdg = std::abs(particle.pdgId());
                const bool is_gen_top = abs_pdg == topPdgId;
                const bool is_gen_boson = bosons.count(abs_pdg) != 0;
                const bool is_lepton = chargedLeptons.count(abs_pdg) || neutralLeptons.count(abs_pdg);
                if((is_gen_top && saveGenTopInfo) || (is_gen_boson && saveGenBosonInfo) || is_lepton) {
                    mothers_to_store.push_back(findStoredMother(particle));
                    particles_to_store.push_back(&particle);
                }
            }
        }

        if(saveGenTopInfo) {
            edm::Handle<TtGenEvent> topGenEvent;
            event.getByToken(topGenEvent_token, topGenEvent);
            if(topGenEvent.isValid()) {
                analysis::GenEventType genEventType = analysis::GenEventType::Other;
                if(topGenEvent->isFullHadronic())
                    genEventType = analysis::GenEventType::TTbar_Hadronic;
                else if(topGenEvent->isSemiLeptonic())
                    genEventType = analysis::GenEventType::TTbar_SemiLeptonic;
                else if(topGenEvent->isFullLeptonic())
                    genEventType = analysis::GenEventType::TTbar_Leptonic;
                (*expressTuple)().genEventType = static_cast<int>(genEventType);
                ++genEventTypeCountMap[genEventType];

                auto top = topGenEvent->top();
                (*expressTuple)().gen_top_pt = top ? top->pt() : ntuple::DefaultFillValue<Float_t>();
                auto top_bar = topGenEvent->topBar();
                (*expressTuple)().gen_topBar_pt = top_bar ? top_bar->pt() : ntuple::DefaultFillValue<Float_t>();
            }
        }

        auto returnIndex = [&](const reco::GenParticle* particle_ptr) {
    		int particle_index = -1;
    		if(particle_ptr != nullptr){
                particle_index = static_cast<int>(particle_ptr - genParticles->data());
    		    if(particle_index > static_cast<int>(genParticles->size()) || particle_index < 0) {
	                  throw analysis::exception("Particle index = %1% for particle with pdgId = %2% exceeds"
                                                " the size of gen particles collection = %3%.")
                                                % particle_index % particle_ptr->pdgId() % genParticles->size();
                }
    		}
    		return particle_index;
    	};

    	for(size_t n = 0; n < particles_to_store.size(); ++n) {
            const reco::GenParticle* particle = particles_to_store.at(n);
            const reco::GenParticle* connected_mother = mothers_to_store.at(n);
            int index = returnIndex(particle);
            (*expressTuple)().genParticles_index.push_back(index);
            (*expressTuple)().genParticles_pdg.push_back(particle->pdgId());
            if(connected_mother) {
                (*expressTuple)().genParticles_rel_pIndex.push_back(index);
                const int mother_index = returnIndex(connected_mother);
                (*expressTuple)().genParticles_rel_mIndex.push_back(mother_index);
            }
        }
    }

private:
    const clock::time_point start;
    const bool isMC, saveGenTopInfo, saveGenBosonInfo;

    edm::EDGetTokenT<LHEEventProduct> lheEventProduct_token;
    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<TtGenEvent> topGenEvent_token;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticles_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;

    std::shared_ptr<ntuple::SummaryTuple> summaryTuple;
    std::shared_ptr<ntuple::ExpressTuple> expressTuple;
    GenEventCountMap genEventCountMap;
    GenEventTypeCountMap genEventTypeCountMap;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SummaryProducer);
