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
#include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"

#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Analysis/include/SummaryTuple.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/TriggerResults.h"
#include "h-tautau/Production/interface/GenTruthTools.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Production/interface/TriggerFileDescriptor.h"
#include "h-tautau/Production/interface/TriggerFileConfigEntryReader.h"
#include "h-tautau/Production/interface/TriggerTools.h"

class SummaryProducer : public edm::EDAnalyzer {
public:
    using clock = std::chrono::system_clock;
    using GenId = ntuple::GenId;
    using GenEventCountMap = ntuple::GenEventCountMap;
    using GenEventTypeCountMap = ntuple::GenEventTypeCountMap;
    using Channel = analysis::Channel;
    using TriggerDescriptorCollection = analysis::TriggerDescriptorCollection;

    SummaryProducer(const edm::ParameterSet& cfg) :
        start(clock::now()),
        isMC(cfg.getParameter<bool>("isMC")),
        saveGenTopInfo(cfg.getParameter<bool>("saveGenTopInfo")),
        lheEventProduct_token(mayConsume<LHEEventProduct>(cfg.getParameter<edm::InputTag>("lheEventProduct"))),
        genEvent_token(mayConsume<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("genEvent"))),
        topGenEvent_token(mayConsume<TtGenEvent>(cfg.getParameter<edm::InputTag>("topGenEvent"))),
        puInfo_token(mayConsume<std::vector<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("puInfo"))),
        taus_token(mayConsume<std::vector<pat::Tau>>(cfg.getParameter<edm::InputTag>("taus"))),
        summaryTuple("summary", &edm::Service<TFileService>()->file(), false)
    {
        summaryTuple().numberOfProcessedEvents = 0;
        if(isMC)
            expressTuple = std::shared_ptr<ntuple::ExpressTuple>(
                    new ntuple::ExpressTuple("all_events", &edm::Service<TFileService>()->file(), false));

        analysis::ConfigReader config_reader;

        analysis::TriggerFileDescriptorCollection trigger_file_descriptors;
        analysis::TriggerFileConfigEntryReader trigger_entry_reader(trigger_file_descriptors);
        config_reader.AddEntryReader("PATTERN", trigger_entry_reader, true);


        const auto& triggerCfg = cfg.getParameter<std::string>("triggerCfg");
        config_reader.ReadConfig(triggerCfg);

        std::map<Channel, TriggerDescriptorCollection> triggerDescriptors;
        const auto& triggerSetup = cfg.getParameterSetVector("triggerSetup");
        for(const auto& channelSetup : triggerSetup) {
            const std::string channel_name = channelSetup.getParameter<std::string>("channel");
            const Channel channel = analysis::Parse<Channel>(channel_name);

            triggerDescriptors[channel] = analysis::TriggerTools::CreateTriggerDescriptors(trigger_file_descriptors,channel);
        }

        for(const auto& channel_desc : triggerDescriptors) {
            const Channel channel = channel_desc.first;
            const int channel_id = static_cast<int>(channel);
            const TriggerDescriptorCollection& descs = channel_desc.second;
            for(size_t n = 0; n < descs.Size(); ++n) {
                summaryTuple().triggers_channel.push_back(channel_id);
                summaryTuple().triggers_pattern.push_back(descs.GetTriggerDescriptor(n).pattern);
            }
        }
    }

private:
    virtual void analyze(const edm::Event& event, const edm::EventSetup&) override
    {
        summaryTuple().numberOfProcessedEvents++;

        if(!tauId_names.size()) {
            edm::Handle<std::vector<pat::Tau>> taus;
            event.getByToken(taus_token, taus);
            for(const auto& tau : *taus) {
                for(const auto& id : tau.tauIDs()) {
                    tauId_names.insert(id.first);
                }
            }
        }

        if(!isMC)
            return;

        edm::Handle<GenEventInfoProduct> genEvent;
        event.getByToken(genEvent_token, genEvent);

        edm::Handle<std::vector<PileupSummaryInfo>> puInfo;
        event.getByToken(puInfo_token, puInfo);
        (*expressTuple)().npu = analysis::gen_truth::GetNumberOfPileUpInteractions(puInfo);
        (*expressTuple)().genEventWeight = genEvent->weight();
        (*expressTuple)().gen_top_pt = ntuple::DefaultFillValue<Float_t>();
        (*expressTuple)().gen_topBar_pt = ntuple::DefaultFillValue<Float_t>();
        (*expressTuple)().lhe_H_m = ntuple::DefaultFillValue<Float_t>();

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
            }
        }

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
        for(const auto& count_entry : genEventTypeCountMap) {
            summaryTuple().genEventType.push_back(static_cast<int>(count_entry.first));
            summaryTuple().genEventType_n_events.push_back(count_entry.second);
        }
        const auto stop = clock::now();
        summaryTuple().exeTime = std::chrono::duration_cast<std::chrono::seconds>(stop - start).count();
        summaryTuple.Fill();
        summaryTuple.Write();
        if(expressTuple)
            expressTuple->Write();
    }

private:
    const clock::time_point start;
    const bool isMC, saveGenTopInfo;

    edm::EDGetTokenT<LHEEventProduct> lheEventProduct_token;
    edm::EDGetTokenT<GenEventInfoProduct> genEvent_token;
    edm::EDGetTokenT<TtGenEvent> topGenEvent_token;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_token;
    edm::EDGetTokenT<std::vector<pat::Tau>> taus_token;

    ntuple::SummaryTuple summaryTuple;
    std::shared_ptr<ntuple::ExpressTuple> expressTuple;
    std::unordered_set<std::string> tauId_names;
    GenEventCountMap genEventCountMap;
    GenEventTypeCountMap genEventTypeCountMap;
};

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SummaryProducer);
