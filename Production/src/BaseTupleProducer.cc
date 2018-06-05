#include "TROOT.h"
#include "h-tautau/Production/interface/BaseTupleProducer.h"
#include "h-tautau/McCorrections/include/TauUncertainties.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "../interface/GenTruthTools.h"
#include "h-tautau/Analysis/include/MetFilters.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/Analysis/include/EventInfo.h"


namespace {
bool EnableThreadSafety() { ROOT::EnableThreadSafety(); return true; }
}

const bool BaseTupleProducer::enableThreadSafety = EnableThreadSafety();

BaseTupleProducer::BaseTupleProducer(const edm::ParameterSet& iConfig, analysis::Channel _channel) :
    treeName(ToString(_channel)),
    anaData(&edm::Service<TFileService>()->file(), treeName + "_stat"),
    electronsMiniAOD_token(mayConsume<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronSrc"))),
    eleTightIdMap_token(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
    eleMediumIdMap_token(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
    tausMiniAOD_token(mayConsume<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
    muonsMiniAOD_token(mayConsume<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    vtxMiniAOD_token(mayConsume<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
    pfMETAOD_token(mayConsume<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
    jetsMiniAOD_token(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
    fatJetsMiniAOD_token(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("fatJetSrc"))),
    PUInfo_token(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
    lheEventProduct_token(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
    genWeights_token(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
    topGenEvent_token(mayConsume<TtGenEvent>(iConfig.getParameter<edm::InputTag>("topGenEvent"))),
    genParticles_token(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    genJets_token(mayConsume<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    productionMode(analysis::EnumNameMap<ProductionMode>::GetDefault().Parse(
                       iConfig.getParameter<std::string>("productionMode"))),
    period(analysis::EnumNameMap<analysis::Period>::GetDefault().Parse(
                       iConfig.getParameter<std::string>("period"))),
    isMC(iConfig.getParameter<bool>("isMC")),
    applyTriggerMatch(iConfig.getParameter<bool>("applyTriggerMatch")),
    runSVfit(iConfig.getParameter<bool>("runSVfit")),
    runKinFit(iConfig.getParameter<bool>("runKinFit")),
    applyRecoilCorr(iConfig.getParameter<bool>("applyRecoilCorr")),
    nJetsRecoilCorr(iConfig.getParameter<int>("nJetsRecoilCorr")),
    saveGenTopInfo(iConfig.getParameter<bool>("saveGenTopInfo")),
    saveGenBosonInfo(iConfig.getParameter<bool>("saveGenBosonInfo")),
    saveGenJetInfo(iConfig.getParameter<bool>("saveGenJetInfo")),
    eventTuple_ptr(ntuple::CreateEventTuple(ToString(_channel),&edm::Service<TFileService>()->file(),false,ntuple::TreeState::Full)),
    eventTuple(*eventTuple_ptr),
    triggerTools(mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "SIM")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "PAT")),
                 consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")),
                 consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")),
                 mayConsume<std::vector<l1extra::L1JetParticle>>(
                                                    iConfig.getParameter<edm::InputTag>("l1JetParticleProduct")),
                 iConfig.getParameter<std::string>("triggerCfg"),
                 _channel)
{
    root_ext::HistogramFactory<TH1D>::LoadConfig(
            edm::FileInPath("h-tautau/Production/data/histograms.cfg").fullPath());
    const std::vector<std::string> energyScaleStrings = iConfig.getParameter<std::vector<std::string>>("energyScales");
    for(const auto& scaleString : energyScaleStrings) {
        const auto es = analysis::Parse<analysis::EventEnergyScale>(scaleString);
        eventEnergyScales.push_back(es);
    }

    if(period == analysis::Period::Run2016){
        badPFMuonFilter_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"));
        badChCandidateFilter_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("badChCandidateFilter"));
    }

    m_rho_token = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));

    if(runSVfit)
//        svfitProducer = std::shared_ptr<analysis::sv_fit::FitProducer>(new analysis::sv_fit::FitProducer(
//            edm::FileInPath("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root").fullPath()));
        svfitProducer = std::shared_ptr<analysis::sv_fit::FitProducer>(new analysis::sv_fit::FitProducer());
    if(runKinFit)
        kinfitProducer = std::shared_ptr<analysis::kin_fit::FitProducer>(new analysis::kin_fit::FitProducer());

    if(applyRecoilCorr)
        recoilPFMetCorrector = std::shared_ptr<RecoilCorrector>(new RecoilCorrector(
            edm::FileInPath("HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root").fullPath()));
}

void BaseTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    InitializeAODCollections(iEvent, iSetup);
    primaryVertex = vertices->ptrAt(0);
    for(auto energyScale : eventEnergyScales) {

        std::ostringstream ss_energyScales;
        ss_energyScales << "events_" << energyScale;
        const std::string energyScales = ss_energyScales.str();
        InitializeCandidateCollections(energyScale);
        try {
            Cutter cut(&GetAnaData().Selection(energyScales));
            cut(true,"events");
            ProcessEvent(cut);
        } catch(cuts::cut_failed&){}

        GetAnaData().Selection(energyScales).fill_selection();
    }
}

void BaseTupleProducer::endJob()
{
    eventTuple.Write();
}

void BaseTupleProducer::InitializeAODCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edmEvent = &iEvent;
    eventId = iEvent.id();
    triggerTools.Initialize(iEvent);
    iEvent.getByToken(electronsMiniAOD_token, pat_electrons);
    iEvent.getByToken(eleTightIdMap_token, tight_id_decisions);
    iEvent.getByToken(eleMediumIdMap_token, medium_id_decisions);
    iEvent.getByToken(tausMiniAOD_token, pat_taus);
    iEvent.getByToken(muonsMiniAOD_token, pat_muons);
    iEvent.getByToken(vtxMiniAOD_token, vertices);
    iEvent.getByToken(pfMETAOD_token, pfMETs);
    iEvent.getByToken(jetsMiniAOD_token, pat_jets);
    iEvent.getByToken(fatJetsMiniAOD_token, pat_fatJets);
    iEvent.getByToken(PUInfo_token, PUInfo);
    if(isMC) {
        iEvent.getByToken(genWeights_token, genEvt);
        iEvent.getByToken(genParticles_token, genParticles);
        iEvent.getByToken(lheEventProduct_token, lheEventProduct);
        iEvent.getByToken(genJets_token, genJets);
        if(saveGenTopInfo)
            iEvent.getByToken(topGenEvent_token, topGenEvent);
    }
    iEvent.getByToken(m_rho_token, rho);

    iSetup.get<JetCorrectionsRecord>().get("AK4PFchs", jetCorParColl);
    jecUnc = std::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty((*jetCorParColl)["Uncertainty"]));

    resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
}

void BaseTupleProducer::InitializeCandidateCollections(analysis::EventEnergyScale energyScale)
{
    using analysis::EventEnergyScale;
    using METUncertainty = pat::MET::METUncertainty;

    static const std::map<EventEnergyScale, int> tauEnergyScales = {
        { EventEnergyScale::TauUp, +1 }, { EventEnergyScale::TauDown, -1 }
    };
    static const std::map<EventEnergyScale, int> jetEnergyScales = {
        { EventEnergyScale::JetUp, +1 }, { EventEnergyScale::JetDown, -1 }
    };
    static const std::map<EventEnergyScale, METUncertainty> metUncertantyMap = {
        { EventEnergyScale::JetUp, METUncertainty::JetEnUp },
        { EventEnergyScale::JetDown, METUncertainty::JetEnDown }
    };

    static const std::map<analysis::Period, std::map<int, double>> tau_correction_factor = {
        { analysis::Period::Run2016, { {0, analysis::uncertainties::tau_2016::sf_1prong},
                                       {1, analysis::uncertainties::tau_2016::sf_1prongPi0},
                                       {10, analysis::uncertainties::tau_2016::sf_3prong} }},
        { analysis::Period::Run2017, { {0, analysis::uncertainties::tau_2017::sf_1prong},
                                       {1, analysis::uncertainties::tau_2017::sf_1prongPi0},
                                       {10, analysis::uncertainties::tau_2017::sf_3prong} } }
    };

    static const std::map<analysis::Period, double> tau_energyUncertainty = {
        { analysis::Period::Run2016, analysis::uncertainties::tau_2016::energyUncertainty},
        { analysis::Period::Run2017, analysis::uncertainties::tau_2017::energyUncertainty}
    };

    eventEnergyScale = energyScale;

    electrons.clear();
    for(size_t n = 0; n < pat_electrons->size(); ++n) {
        const edm::Ptr<pat::Electron> ele_ptr(pat_electrons, n);
        electrons.push_back(ElectronCandidate(ele_ptr, Isolation(*ele_ptr)));
    }

    muons.clear();
    for(const auto& muon : *pat_muons)
        muons.push_back(MuonCandidate(muon, Isolation(muon)));


    taus.clear();
    for(const auto& tau : *pat_taus) {
        TauCandidate tauCandidate(tau, Isolation(tau));
        if(isMC) {
            const analysis::gen_truth::MatchResult result =
                        analysis::gen_truth::LeptonGenMatch(tauCandidate.GetMomentum(), *genParticles);
            if(result.first == analysis::GenMatch::Tau){
                double corr_factor = 1;
                double sf = 1;
                if(tau_correction_factor.at(period).count(tau.decayMode())){
                    corr_factor = tau_correction_factor.at(period).at(tau.decayMode());
                    sf = corr_factor;
                    if(tauEnergyScales.count(energyScale)) {
                        const int sign = tauEnergyScales.at(energyScale);
                        sf = corr_factor + sign * tau_energyUncertainty.at(period);
                    }
                }
                const auto shiftedMomentum = tau.p4() * sf;
                tauCandidate.SetMomentum(shiftedMomentum);
            }
        }

        taus.push_back(tauCandidate);
    }

    jets.clear();
    for(const auto& jet : *pat_jets) {
        JetCandidate jetCandidate(jet);
        if(jetEnergyScales.count(energyScale)) {
            jecUnc->setJetEta(jet.eta());
            jecUnc->setJetPt(jet.pt()); // here you must use the CORRECTED jet pt
            const double unc = jecUnc->getUncertainty(true);
            const int sign = jetEnergyScales.at(energyScale);
            const double sf = (1.0 + (sign * unc));
            const auto shiftedMomentum = jet.p4() * sf;
            jetCandidate.SetMomentum(shiftedMomentum);
        }
        jets.push_back(jetCandidate);
    }

    fatJets.clear();
    for(const auto& jet : * pat_fatJets)
        fatJets.push_back(JetCandidate(jet));

    met = std::shared_ptr<MET>(new MET((*pfMETs)[0], (*pfMETs)[0].getSignificanceMatrix()));
    if(metUncertantyMap.count(energyScale)) {
        const auto shiftedMomentum = (*met)->shiftedP4(metUncertantyMap.at(energyScale));
        met->SetMomentum(shiftedMomentum);
    }
}

const double BaseTupleProducer::Isolation(const pat::Electron& electron)
{
    const double sum_neutral = electron.pfIsolationVariables().sumNeutralHadronEt
                             + electron.pfIsolationVariables().sumPhotonEt
                             - 0.5 * electron.pfIsolationVariables().sumPUPt;
    const double abs_iso = electron.pfIsolationVariables().sumChargedHadronPt + std::max(sum_neutral, 0.0);
    return abs_iso / electron.pt();
}

const double BaseTupleProducer::Isolation(const pat::Muon& muon)
{
    const double sum_neutral = muon.pfIsolationR04().sumNeutralHadronEt
                             + muon.pfIsolationR04().sumPhotonEt
                             - 0.5 * muon.pfIsolationR04().sumPUPt;
    const double abs_iso = muon.pfIsolationR04().sumChargedHadronPt + std::max(sum_neutral, 0.0);
    return abs_iso / muon.pt();
}

const double BaseTupleProducer::Isolation(const pat::Tau& tau)
{
    static const std::map<analysis::Period, analysis::TauIdDiscriminator> discriminators = {
        { analysis::Period::Run2016, analysis::TauIdDiscriminator::byIsolationMVArun2v1DBoldDMwLT2016 },
        { analysis::Period::Run2017, analysis::TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017 }
    };
    const auto& desc = analysis::tau_id::GetTauIdDescriptors().at(discriminators.at(period));
    return tau.tauID(desc.ToStringRaw());
}

analysis::TauIdResults BaseTupleProducer::CreateTauIdResults(const pat::Tau& tau, analysis::Period period)
{
    using analysis::Period;
    using analysis::TauIdDiscriminator;
    using analysis::TauIdResults;
    static const std::map<Period, std::set<TauIdDiscriminator>> ignore_list {
        { Period::Run2016, { TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMwLT2017,
                             TauIdDiscriminator::byIsolationMVArun2017v2DBoldDMdR0p3wLT2017 } },
        { Period::Run2017, {} }
    };
    analysis::TauIdResults results;
    const auto& descs = analysis::TauIdResults::GetResultDescriptors();
    const auto& ignore = ignore_list.at(period);
    for(size_t n = 0; n < descs.size(); ++n) {
        if(!ignore.count(descs.at(n).discriminator))
            results.SetResult(n, tau.tauID(descs.at(n).ToString()) > .5f);
    }
    return results;
}

//https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017#Preliminary_Recommendations_for
//recommended for 2017
bool BaseTupleProducer::PassPFTightId(const pat::Jet& pat_jet, analysis::Period period)
{
    const pat::Jet& patJet = pat_jet.correctedJet("Uncorrected");
    const double abs_eta = std::abs(patJet.eta());
    if(period == analysis::Period::Run2016)
    {
        if(abs_eta <= 2.7 && (
                         patJet.neutralHadronEnergyFraction() >= 0.9 ||
                         patJet.neutralEmEnergyFraction() >= 0.9 ||
                         patJet.nConstituents() <= 1)) return false;
        if(abs_eta <= 2.4 && (
                          patJet.chargedHadronEnergyFraction() <= 0 ||
                          patJet.chargedMultiplicity() <= 0 ||
                          patJet.chargedEmEnergyFraction() >= 0.99)) return false;

        if(abs_eta > 2.7 && abs_eta <= 3.0 && (
                                           patJet.neutralEmEnergyFraction() <= 0.01 ||
                                           patJet.neutralEmEnergyFraction() >= 0.99 ||
                                           patJet.neutralMultiplicity() <= 2)) return false;

        if(abs_eta > 3.0 && (
                         patJet.neutralEmEnergyFraction() >= 0.9 ||
                         patJet.neutralMultiplicity() <= 10)) return false;
    }

    if(period == analysis::Period::Run2017)
    {
        if(abs_eta <= 2.7 && (
                         patJet.neutralHadronEnergyFraction() >= 0.9 ||
                         patJet.neutralEmEnergyFraction() >= 0.9 ||
                         patJet.nConstituents() <= 1)) return false;
        if(abs_eta <= 2.4 && (
                          patJet.chargedHadronEnergyFraction() <= 0 ||
                          patJet.chargedMultiplicity() <= 0 )) return false;

        if(abs_eta > 2.7 && abs_eta <= 3.0 && (
                                           patJet.neutralEmEnergyFraction() <= 0.02 ||
                                           patJet.neutralEmEnergyFraction() >= 0.99 ||
                                           patJet.neutralMultiplicity() <= 2)) return false;

        if(abs_eta > 3.0 && (
                         patJet.neutralEmEnergyFraction() >= 0.9 ||
                         patJet.neutralHadronEnergyFraction() <= 0.02 ||
                         patJet.neutralMultiplicity() <= 10)) return false;

    }
    return true;
}

void BaseTupleProducer::FillLheInfo(bool haveReference)
{
    if(haveReference || !lheEventProduct.isValid()) {
        eventTuple().lhe_n_partons = ntuple::DefaultFillValue<UInt_t>();
        eventTuple().lhe_n_c_partons = ntuple::DefaultFillValue<UInt_t>();
        eventTuple().lhe_n_b_partons = ntuple::DefaultFillValue<UInt_t>();
        eventTuple().lhe_HT = ntuple::DefaultFillValue<Float_t>();
        eventTuple().lhe_H_m = ntuple::DefaultFillValue<Float_t>();
        eventTuple().lhe_hh_m = ntuple::DefaultFillValue<Float_t>();
        eventTuple().lhe_hh_cosTheta = ntuple::DefaultFillValue<Float_t>();
        return;
    }

    const auto lheSummary = analysis::gen_truth::ExtractLheSummary(*lheEventProduct);
    eventTuple().lhe_n_partons = lheSummary.n_partons;
    eventTuple().lhe_n_c_partons = lheSummary.n_c_partons;
    eventTuple().lhe_n_b_partons = lheSummary.n_b_partons;
    eventTuple().lhe_HT = lheSummary.HT;
    eventTuple().lhe_H_m = lheSummary.m_H;
    eventTuple().lhe_hh_m = lheSummary.m_hh;
    eventTuple().lhe_hh_cosTheta = lheSummary.cosTheta_hh;
}

void BaseTupleProducer::ApplyRecoilCorrection(const std::vector<JetCandidate>& jets)
{
    analysis::LorentzVectorXYZ total_p4, vis_p4;
    for(const auto& particle : *genParticles) {
        if( (particle.fromHardProcessFinalState() && (particle.isMuon() || particle.isElectron()|| reco::isNeutrino(particle)))
              ||particle.isDirectHardProcessTauDecayProductFinalState())
            total_p4 += particle.p4();
        if( (particle.fromHardProcessFinalState() && (particle.isMuon() || particle.isElectron()))
	            || (particle.isDirectHardProcessTauDecayProductFinalState() && !reco::isNeutrino(particle)) )
            vis_p4 += particle.p4();
    }


    const double genPx = total_p4.px();
    const double genPy = total_p4.py();

    const double visPx = vis_p4.px();
    const double visPy = vis_p4.py();

    float pfmetcorr_ex, pfmetcorr_ey;

    const double pfmet_ex = met->GetMomentum().px();
    const double pfmet_ey = met->GetMomentum().py();

    int njets= nJetsRecoilCorr;
    for(const auto& jet : jets)
    {
        if(jet.GetMomentum().pt() > 30) njets++;
    }
    recoilPFMetCorrector->CorrectByMeanResolution(
	    pfmet_ex, // uncorrected type I pf met px (float)
	    pfmet_ey, // uncorrected type I pf met py (float)
	    genPx, // generator Z/W/Higgs px (float)
	    genPy, // generator Z/W/Higgs py (float)
	    visPx, // generator visible Z/W/Higgs px (float)
	    visPy, // generator visible Z/W/Higgs py (float)
	    njets,  // number of jets (hadronic jet multiplicity) (int)
	    pfmetcorr_ex, // corrected type I pf met px (float)
	    pfmetcorr_ey  // corrected type I pf met py (float)
	    );

   const TVector2 met_vect(pfmetcorr_ex,pfmetcorr_ey);
   met->SetMomentum(analysis::LorentzVectorM(met_vect.Mod(),0,met_vect.Phi(),0));
}

void BaseTupleProducer::FillGenParticleInfo()
{
    using analysis::GenEventType;
    static constexpr int electronPdgId = 11, muonPdgId = 13, tauPdgId = 15;
    static const std::set<int> bosons = { 23, 24, 25, 35 };

    std::vector<const reco::GenParticle*> particles_to_store;

    std::map<int, size_t> particle_counts;
    for(const auto& particle : *genParticles) {
        const auto& flag = particle.statusFlags();
        if(!flag.isPrompt() || !flag.isLastCopy()) continue;
        const int abs_pdg = std::abs(particle.pdgId());
        ++particle_counts[abs_pdg];
        if(saveGenBosonInfo && bosons.count(abs_pdg))
            particles_to_store.push_back(&particle);
    }

    eventTuple().genParticles_nPromptElectrons = particle_counts[electronPdgId];
    eventTuple().genParticles_nPromptMuons = particle_counts[muonPdgId];
    eventTuple().genParticles_nPromptTaus = particle_counts[tauPdgId];

    if(saveGenTopInfo) {
        GenEventType genEventType = GenEventType::Other;
        if(topGenEvent->isFullHadronic())
                    genEventType = GenEventType::TTbar_Hadronic;
        else if(topGenEvent->isSemiLeptonic())
            genEventType = GenEventType::TTbar_SemiLeptonic;
        else if(topGenEvent->isFullLeptonic())
            genEventType = GenEventType::TTbar_Leptonic;
        eventTuple().genEventType = static_cast<int>(genEventType);

        auto top = topGenEvent->top();
        if(top)
            particles_to_store.push_back(top);
        auto top_bar = topGenEvent->topBar();
        if(top_bar)
            particles_to_store.push_back(top_bar);
    }

    for(auto particle : particles_to_store) {
        eventTuple().genParticles_pdg.push_back(particle->pdgId());
        eventTuple().genParticles_p4.push_back(ntuple::LorentzVectorM(particle->p4()));
    }

}

void BaseTupleProducer::FillGenJetInfo()
{
    static constexpr int b_flavour = 5, c_flavour = 4;
    static constexpr double pt_cut = 5;
    eventTuple().genJets_nTotal = genJets->size();

    std::map<int, size_t> pf_counts, hf_counts;

    for(const JetCandidate& jet : jets) {
        ++pf_counts[std::abs(jet->partonFlavour())];
        ++hf_counts[std::abs(jet->hadronFlavour())];
    }
    eventTuple().jets_nTotal_partonFlavour_b = pf_counts[b_flavour];
    eventTuple().jets_nTotal_partonFlavour_c = pf_counts[c_flavour];
    eventTuple().jets_nTotal_hadronFlavour_b = hf_counts[b_flavour];
    eventTuple().jets_nTotal_hadronFlavour_c = hf_counts[c_flavour];

    if(!saveGenJetInfo) return;

    for(const reco::GenJet& gen_jet : *genJets) {
        if(gen_jet.pt() <= pt_cut) continue;
        eventTuple().genJets_p4.push_back(ntuple::LorentzVectorE(gen_jet.p4()));

        const auto findRecoJetFlavours = [&]() -> std::pair<int, int> {
            for(const JetCandidate& reco_jet : jets) {
                if(reco_jet->genJet() == &gen_jet)
                    return std::make_pair(reco_jet->partonFlavour(), reco_jet->hadronFlavour());
            }
            return std::make_pair(ntuple::DefaultFillValue<int>(), ntuple::DefaultFillValue<int>());
        };

        const auto flavours = findRecoJetFlavours();
        eventTuple().genJets_partonFlavour.push_back(flavours.first);
        eventTuple().genJets_hadronFlavour.push_back(flavours.second);
    }
}

void BaseTupleProducer::FillLegGenMatch(size_t leg_id, const analysis::LorentzVectorXYZ& p4)
{
    using namespace analysis;
    static constexpr int default_int_value = ntuple::DefaultFillValue<Int_t>();

    auto& gen_match = leg_id == 1 ? eventTuple().gen_match_1 : eventTuple().gen_match_2;
    auto& gen_p4 = leg_id == 1 ? eventTuple().gen_p4_1 : eventTuple().gen_p4_2;

    if(isMC) {
        const auto match = gen_truth::LeptonGenMatch(p4, *genParticles);
        gen_match = static_cast<int>(match.first);
        const auto metched_p4 = match.second ? match.second->p4() : LorentzVectorXYZ();
        gen_p4 = ntuple::LorentzVectorM(metched_p4);
    } else {
        gen_match = default_int_value;
        gen_p4 = ntuple::LorentzVectorM();
    }
}

void BaseTupleProducer::FillMetFilters(analysis::Period period)
{
    using MetFilters = ntuple::MetFilters;
    using Filter = MetFilters::Filter;

    MetFilters filters;
    const auto setResult = [&](Filter filter, const std::string& name) {
        bool result;
        if(!triggerTools.TryGetAnyTriggerResult(name, result))
            result = true;
        filters.SetResult(filter, result);
    };

    setResult(Filter::PrimaryVertex, "Flag_goodVertices");
    setResult(Filter::BeamHalo, "Flag_globalTightHalo2016Filter");
    setResult(Filter::HBHE_noise, "Flag_HBHENoiseFilter");
    setResult(Filter::HBHEiso_noise, "Flag_HBHENoiseIsoFilter");
    setResult(Filter::ECAL_TP, "Flag_EcalDeadCellTriggerPrimitiveFilter");
    setResult(Filter::ee_badSC_noise, "Flag_eeBadScFilter");

    if(period == analysis::Period::Run2016){
        edm::Handle<bool> badPFMuon;
        edmEvent->getByToken(badPFMuonFilter_token, badPFMuon);
        filters.SetResult(Filter::badMuon, *badPFMuon);

        edm::Handle<bool> badChCandidate;
        edmEvent->getByToken(badChCandidateFilter_token, badChCandidate);
        filters.SetResult(Filter::badChargedHadron,*badChCandidate);
    }

    if(period == analysis::Period::Run2017){
        setResult(Filter::badMuon, "Flag_BadPFMuonFilter");
        setResult(Filter::badChargedHadron, "Flag_BadChargedCandidateFilter");
        setResult(Filter::ecalBadCalib, "Flag_ecalBadCalibFilter");
    }


    eventTuple().metFilters = filters.FilterResults();
}

void BaseTupleProducer::ApplyBaseSelection(analysis::SelectionResultsBase& selection,
                        const std::vector<LorentzVector>& signalLeptonMomentums)
{
    static constexpr bool RunKinfitForAllPairs = false;
    using namespace cuts::btag_2016;

    selection.jets = CollectJets(signalLeptonMomentums);
    const size_t n_jets = selection.jets.size();
    if(applyRecoilCorr)
        ApplyRecoilCorrection(selection.jets);
    if(!runKinFit || n_jets < 2) return;

    const auto runKinfit = [&](size_t first, size_t second) {
        const ntuple::JetPair pair(first, second);
        //jet resolution
        JME::JetParameters parameters_1;
        parameters_1.setJetPt(selection.jets.at(pair.first).GetMomentum().pt());
        parameters_1.setJetEta(selection.jets.at(pair.first).GetMomentum().eta());
        parameters_1.setRho(*rho);
        float resolution_1 = resolution.getResolution(parameters_1);
        float energy_resolution_1 = resolution_1 * selection.jets.at(pair.first).GetMomentum().E();
        JME::JetParameters parameters_2;
        parameters_2.setJetPt(selection.jets.at(pair.second).GetMomentum().pt());
        parameters_2.setJetEta(selection.jets.at(pair.second).GetMomentum().eta());
        parameters_2.setRho(*rho);
        float resolution_2 = resolution.getResolution(parameters_2);
        float energy_resolution_2 = resolution_2 * selection.jets.at(pair.second).GetMomentum().E();

        const size_t pair_index = ntuple::CombinationPairToIndex(pair, n_jets);
        const auto& result = kinfitProducer->Fit(signalLeptonMomentums.at(0), signalLeptonMomentums.at(1),
                                                 selection.jets.at(pair.first).GetMomentum(),
                                                 selection.jets.at(pair.second).GetMomentum(), *met,energy_resolution_1,
                                                 energy_resolution_2);
        selection.kinfitResults[pair_index] = result;
    };

    if(RunKinfitForAllPairs) {
        for(size_t n = 0; n < n_jets; ++n) {
            for(size_t k = 0; k < n_jets; ++k) {
                if(k == n) continue;
                runKinfit(n, k);
            }
        }
    } else {
        runKinfit(0, 1);

        std::vector<analysis::jet_ordering::JetInfo<LorentzVector>> jet_info_vector;
        for(size_t n = 0; n < selection.jets.size(); ++n) {
            const JetCandidate& jet = selection.jets.at(n);
            jet_info_vector.emplace_back(jet.GetMomentum(),n,jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        }

        auto jets_ordered = analysis::jet_ordering::OrderJets(jet_info_vector,false,cuts::btag_2017::pt,
                                                              cuts::btag_2017::eta);
        if((jets_ordered.at(0).index) != 0 && (jets_ordered.at(1).index) != 1){
            runKinfit(jets_ordered.at(0).index, jets_ordered.at(1).index);
        }
    }
}

std::vector<BaseTupleProducer::ElectronCandidate> BaseTupleProducer::CollectVetoElectrons(
        const std::vector<const ElectronCandidate*>& signalElectrons)
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectVetoElectron, this, _1, _2, signalElectrons);
    return CollectObjects("vetoElectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::MuonCandidate> BaseTupleProducer::CollectVetoMuons(
        const std::vector<const MuonCandidate*>& signalMuons)
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectVetoMuon, this, _1, _2, signalMuons);
    return CollectObjects("vetoMuons", base_selector, muons);
}

std::vector<BaseTupleProducer::JetCandidate> BaseTupleProducer::CollectJets(
        const std::vector<LorentzVector>& signalLeptonMomentums)
{
    using namespace cuts::btag_2016;
    using namespace std::placeholders;
    const auto baseSelector = std::bind(&BaseTupleProducer::SelectJet, this, _1, _2, signalLeptonMomentums);

    const auto comparitor = [](const JetCandidate& j1, const JetCandidate& j2) {
        const auto deepcsv1 = j1->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll");
        const auto deepcsv2 = j2->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll");
        const analysis::jet_ordering::JetInfo<LorentzVector> jet_info_1(j1.GetMomentum(),0,deepcsv1);
        const analysis::jet_ordering::JetInfo<LorentzVector> jet_info_2(j2.GetMomentum(),1,deepcsv2);
        return analysis::jet_ordering::CompareJets(jet_info_1,jet_info_2,cuts::btag_2016::pt,cuts::btag_2016::eta);
    };

    return CollectObjects("jets", baseSelector, jets, comparitor);
}

void BaseTupleProducer::SelectVetoElectron(const ElectronCandidate& electron, Cutter& cut,
                                           const std::vector<const ElectronCandidate*>& signalElectrons) const
{
    using namespace cuts::H_tautau_2016::electronVeto;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electron_dxy = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electron_dxy < dxy, "dxy", electron_dxy);
    const double electron_dz = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electron_dz < dz, "dz", electron_dz);
    const bool isMedium  = (*medium_id_decisions)[electron.getPtr()];
    cut(isMedium, "electronMVAMediumID");
    if(productionMode != ProductionMode::hh) {
        cut(electron->passConversionVeto(), "conversionVeto");
    }
    if(period != analysis::Period::Run2017)
        cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
    for(size_t n = 0; n < signalElectrons.size(); ++n) {
        std::ostringstream ss_name;
        ss_name << "isNotSignal_" << n + 1;
        const bool isNotSignal =  &(*electron) != &(*(*signalElectrons.at(n)));
        cut(isNotSignal, ss_name.str());
    }
}

void BaseTupleProducer::SelectVetoMuon(const MuonCandidate& muon, Cutter& cut,
                                       const std::vector<const MuonCandidate*>& signalMuons) const
{
    using namespace cuts::H_tautau_2016::muonVeto;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    cut(muon.GetIsolation() < pfRelIso04, "iso", muon.GetIsolation());

    bool passMuonId =  muon->isMediumMuon();
    if( productionMode == ProductionMode::hh) passMuonId = muon->isLooseMuon();
    cut(passMuonId, "muonID");
    for(size_t n = 0; n < signalMuons.size(); ++n) {
        std::ostringstream ss_name;
        ss_name << "isNotSignal_" << n + 1;
        const bool isNotSignal =  &(*muon) != &(*(*signalMuons.at(n)));
        cut(isNotSignal, ss_name.str());
    }
}

void BaseTupleProducer::SelectJet(const JetCandidate& jet, Cutter& cut,
                                  const std::vector<LorentzVector>& signalLeptonMomentums) const
{
    using namespace cuts::H_tautau_2016::jetID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = jet.GetMomentum();
    cut(p4.Pt() > pt - pt_safety, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < cuts::hh_bbtautau_2017::jetID::eta, "eta", p4.Eta());
    cut(PassPFTightId(*jet,period), "jet_id");
    for(size_t n = 0; n < signalLeptonMomentums.size(); ++n) {
        std::ostringstream cut_name;
        cut_name << "deltaR_lep" << n + 1;
        const double deltaR = ROOT::Math::VectorUtil::DeltaR(p4, signalLeptonMomentums.at(n));
        cut(deltaR > deltaR_signalObjects, cut_name.str(), deltaR);
    }
}


#define GET_LEG(x) (leg_id == 1 ? eventTuple().x##_1 : eventTuple().x##_2)
#define RAW_ID(name, n) GET_LEG(tauId_##name) = fill_tauIds && !ignore.count(#name) ? tau->tauID(#name) : default_value;

void BaseTupleProducer::FillElectronLeg(size_t leg_id, const ElectronCandidate& electron)
{
    GET_LEG(p4) = ntuple::LorentzVectorM(electron.GetMomentum());
    GET_LEG(q) = electron.GetCharge();
    GET_LEG(dxy) = electron->gsfTrack()->dxy(primaryVertex->position());
    GET_LEG(dz) = electron->gsfTrack()->dz(primaryVertex->position());
    GET_LEG(iso) = electron.GetIsolation();
    GET_LEG(decayMode) = -1;
    FillLegGenMatch(leg_id, electron->p4());
}

void BaseTupleProducer::FillMuonLeg(size_t leg_id, const MuonCandidate& muon)
{
    GET_LEG(p4) = ntuple::LorentzVectorM(muon.GetMomentum());
    GET_LEG(q) = muon.GetCharge();
    GET_LEG(dxy) = muon->muonBestTrack()->dxy(primaryVertex->position());
    GET_LEG(dz) = muon->muonBestTrack()->dz(primaryVertex->position());
    GET_LEG(iso) = muon.GetIsolation();
    GET_LEG(decayMode) = -1;
    FillLegGenMatch(leg_id, muon->p4());
}

void BaseTupleProducer::FillTauLeg(size_t leg_id, const TauCandidate& tau, bool fill_tauIds)
{
    static const float default_value = ntuple::DefaultFillValue<float>();
    using analysis::Period;
    static const std::map<Period, std::set<std::string>> ignore_list {
        { Period::Run2016, { "byIsolationMVArun2017v2DBoldDMwLTraw2017",
                             "byIsolationMVArun2017v2DBoldDMdR0p3wLTraw2017" } },
        { Period::Run2017, { "againstElectronMVA6RawNew", "againstElectronMVA6categoryNew" } }
    };
    const auto& ignore = ignore_list.at(period);

    GET_LEG(p4) = ntuple::LorentzVectorM(tau.GetMomentum());
    GET_LEG(q) = tau.GetCharge();
    const auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    GET_LEG(dxy) = packedLeadTauCand->dxy();
    GET_LEG(dz) = packedLeadTauCand->dz();
    GET_LEG(iso) = tau.GetIsolation();
    GET_LEG(decayMode) = tau->decayMode();
    FillLegGenMatch(leg_id, tau->p4());
    GET_LEG(tauId_flags) = CreateTauIdResults(*tau, period).GetResultBits();
    RAW_TAU_IDS(leg_id)
}

#undef GET_LEG
#undef RAW_ID

void BaseTupleProducer::FillEventTuple(const analysis::SelectionResultsBase& selection,
                                       const analysis::SelectionResultsBase* reference)
{
    using namespace analysis;
    using EventPart = ntuple::StorageMode::EventPart;

    ntuple::StorageMode storageMode;

    eventTuple().run  = edmEvent->id().run();
    eventTuple().lumi = edmEvent->id().luminosityBlock();
    eventTuple().evt  = edmEvent->id().event();
    eventTuple().eventEnergyScale = static_cast<int>(eventEnergyScale);
    eventTuple().genEventType = static_cast<int>(GenEventType::Other);
    eventTuple().genEventWeight = isMC ? genEvt->weight() : 1;

    eventTuple().npv = vertices->size();
    eventTuple().npu = gen_truth::GetNumberOfPileUpInteractions(PUInfo);
    eventTuple().rho = *rho;

    // HTT candidate
    eventTuple().SVfit_p4 = selection.svfitResult.momentum;
    eventTuple().SVfit_p4_error = selection.svfitResult.momentum_error;
    eventTuple().SVfit_mt = selection.svfitResult.transverseMass;
    eventTuple().SVfit_mt_error = selection.svfitResult.transverseMass_error;

    // MET
    eventTuple().pfMET_p4 = met->GetMomentum();
    eventTuple().pfMET_cov = met->GetCovMatrix();
    FillMetFilters(period);

    std::set<const pat::Jet*> selected_jets;
    if(!reference || !selection.HaveSameJets(*reference)) {
        for(const JetCandidate& jet : selection.jets) {
            const auto selected_jet = &(*jet);
            selected_jets.insert(selected_jet);
            const LorentzVector& p4 = jet.GetMomentum();
            eventTuple().jets_p4.push_back(ntuple::LorentzVectorE(p4));
            eventTuple().jets_csv.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            eventTuple().jets_deepCsv_BvsAll.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")); //sum of b and bb
            eventTuple().jets_deepCsv_CvsB.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB"));
            eventTuple().jets_deepCsv_CvsL.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL"));
            eventTuple().jets_deepFlavour_b.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probb"));
            eventTuple().jets_deepFlavour_bb.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
            eventTuple().jets_deepFlavour_lepb.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
            eventTuple().jets_deepFlavour_c.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probc"));
            eventTuple().jets_deepFlavour_uds.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
            eventTuple().jets_deepFlavour_g.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probg"));
            eventTuple().jets_rawf.push_back((jet->correctedJet("Uncorrected").pt() ) / p4.Pt());
            eventTuple().jets_pu_id.push_back(jet->userInt("pileupJetId:fullId"));
            eventTuple().jets_partonFlavour.push_back(jet->partonFlavour());
            eventTuple().jets_hadronFlavour.push_back(jet->hadronFlavour());
            // Jet resolution
            JME::JetParameters parameters;
            parameters.setJetPt(jet.GetMomentum().pt());
            parameters.setJetEta(jet.GetMomentum().eta());
            parameters.setRho(*rho);
            float jet_resolution = resolution.getResolution(parameters);
            eventTuple().jets_resolution.push_back(jet_resolution); // percentage

            eventTuple().jets_triggerFilterMatch.push_back(triggerTools.GetJetMatchBits(p4,
                                                           cuts::H_tautau_2016::DeltaR_triggerMatch));
        }
        if(eventEnergyScale == EventEnergyScale::Central){
            for(const auto jet_cand : jets){
                const auto pat_jet = &(*jet_cand);
                if(selected_jets.count(pat_jet)) continue;
                const LorentzVector& other_p4 = jet_cand.GetMomentum();
                eventTuple().other_jets_p4.push_back(ntuple::LorentzVectorE(other_p4));
            }
        }
    } else
        storageMode.SetPresence(EventPart::Jets, false);

    const bool haveReference = reference && selection.eventId == reference->eventId;
    storageMode.SetPresence(EventPart::FatJets, !haveReference);
    storageMode.SetPresence(EventPart::GenInfo, !haveReference);
    eventTuple().storageMode = storageMode.Mode();

    if(!haveReference) {
        for(const JetCandidate& jet : fatJets) {
            const LorentzVector& p4 = jet.GetMomentum();
            eventTuple().fatJets_p4.push_back(ntuple::LorentzVectorE(p4));
            eventTuple().fatJets_csv.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            eventTuple().fatJets_deepCsv_BvsAll.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")); //sum of b and bb
            eventTuple().fatJets_deepCsv_CvsB.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB"));
            eventTuple().fatJets_deepCsv_CvsL.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL"));
            eventTuple().fatJets_deepFlavour_b.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probb"));
            eventTuple().fatJets_deepFlavour_bb.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
            eventTuple().fatJets_deepFlavour_lepb.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
            eventTuple().fatJets_deepFlavour_c.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probc"));
            eventTuple().fatJets_deepFlavour_uds.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
            eventTuple().fatJets_deepFlavour_g.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probg"));
            eventTuple().fatJets_m_pruned.push_back(GetUserFloat(jet, "ak8PFJetsCHSPrunedMass"));
            eventTuple().fatJets_m_softDrop.push_back(GetUserFloat(jet, "ak8PFJetsCHSSoftDropMass"));
            eventTuple().fatJets_n_subjettiness_tau1.push_back(GetUserFloat(jet, "NjettinessAK8:tau1"));
            eventTuple().fatJets_n_subjettiness_tau2.push_back(GetUserFloat(jet, "NjettinessAK8:tau2"));
            eventTuple().fatJets_n_subjettiness_tau3.push_back(GetUserFloat(jet, "NjettinessAK8:tau3"));

            if(!jet->hasSubjets("SoftDrop")) continue;
            const size_t parentIndex = eventTuple().fatJets_p4.size() - 1;
            const auto& sub_jets = jet->subjets("SoftDrop");
            for(const auto& sub_jet : sub_jets) {
                eventTuple().subJets_p4.push_back(ntuple::LorentzVectorE(sub_jet->p4()));
                eventTuple().subJets_csv.push_back(
                            sub_jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                eventTuple().subJets_deepCsv_BvsAll.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll")); //sum of b and bb
                eventTuple().subJets_deepCsv_CvsB.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB"));
                eventTuple().subJets_deepCsv_CvsL.push_back(jet->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL"));
                eventTuple().subJets_deepFlavour_b.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probb"));
                eventTuple().subJets_deepFlavour_bb.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probbb"));
                eventTuple().subJets_deepFlavour_lepb.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:problepb"));
                eventTuple().subJets_deepFlavour_c.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probc"));
                eventTuple().subJets_deepFlavour_uds.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probuds"));
                eventTuple().subJets_deepFlavour_g.push_back(jet->bDiscriminator("pfDeepFlavourJetTags:probg"));
                eventTuple().subJets_parentIndex.push_back(parentIndex);
            }
        }
    }

    for(const auto& result : selection.kinfitResults) {
        eventTuple().kinFit_jetPairId.push_back(result.first);
        eventTuple().kinFit_m.push_back(result.second.mass);
        eventTuple().kinFit_chi2.push_back(result.second.chi2);
        eventTuple().kinFit_convergence.push_back(result.second.convergence);
    }

    FillLheInfo(haveReference);
    if(isMC && !haveReference) {
        FillGenParticleInfo();
        FillGenJetInfo();
    }

    eventTuple().extraelec_veto = selection.electronVeto;
    eventTuple().extramuon_veto = selection.muonVeto;

    eventTuple().trigger_match = !applyTriggerMatch || selection.triggerResults.AnyAcceptAndMatch();
    eventTuple().trigger_accepts = selection.triggerResults.GetAcceptBits();
    eventTuple().trigger_matches = selection.triggerResults.GetMatchBits();
}
