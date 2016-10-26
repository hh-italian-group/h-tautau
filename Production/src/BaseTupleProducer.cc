#include "h-tautau/Production/interface/BaseTupleProducer.h"
#include "h-tautau/McCorrections/include/TauUncertainties.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "../interface/GenTruthTools.h"

BaseTupleProducer::BaseTupleProducer(const edm::ParameterSet& iConfig, const std::string& treeName):
    anaData(&edm::Service<TFileService>()->file(), treeName + "_stat"),
    anaDataBeforeCut(&edm::Service<TFileService>()->file(), treeName + "_before_cut"),
    anaDataAfterCut(&edm::Service<TFileService>()->file(), treeName + "_after_cut"),
    anaDataFinalSelection(&edm::Service<TFileService>()->file(), treeName + "_final_selection"),
    electronsMiniAOD_token(mayConsume<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronSrc"))),
    eleTightIdMap_token(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap"))),
    eleMediumIdMap_token(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
    eleCutBasedVetoMap_token(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleCutBasedVeto"))),
    tausMiniAOD_token(mayConsume<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
    muonsMiniAOD_token(mayConsume<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    vtxMiniAOD_token(mayConsume<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
    pfMETAOD_token(mayConsume<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
    jetsMiniAOD_token(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
    fatJetsMiniAOD_token(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("fatJetSrc"))),
    metCovMatrix_token(consumes<MetCovMatrix>(iConfig.getParameter<edm::InputTag>("metCov"))),
    PUInfo_token(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
    lheEventProduct_token(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
    genWeights_token(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
    prunedGen_token(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    productionMode(analysis::EnumNameMap<ProductionMode>::GetDefault().Parse(
                       iConfig.getParameter<std::string>("productionMode"))),
    isMC(iConfig.getParameter<bool>("isMC")),
    applyTriggerMatch(iConfig.getParameter<bool>("applyTriggerMatch")),
    runSVfit(iConfig.getParameter<bool>("runSVfit")),
    runKinFit(iConfig.getParameter<bool>("runKinFit")),
    hltPaths(iConfig.getParameter<std::vector<std::string>>("hltPaths")),
    eventTuple(treeName, &edm::Service<TFileService>()->file(), false),
    triggerTools(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits")),
                 consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")),
                 consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")),
                 mayConsume<std::vector<l1extra::L1JetParticle>>(
                                                    iConfig.getParameter<edm::InputTag>("l1JetParticleProduct"))),
    svfitProducer(edm::FileInPath("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root").fullPath())
{
    root_ext::HistogramFactory<TH1D>::LoadConfig(
            edm::FileInPath("h-tautau/Production/data/histograms.cfg").fullPath());
    const std::vector<std::string> energyScaleStrings = iConfig.getParameter<std::vector<std::string>>("energyScales");
    for(const auto& scaleString : energyScaleStrings) {
        const auto es = analysis::Parse<analysis::EventEnergyScale>(scaleString);
        eventEnergyScales.push_back(es);
    }
}

void BaseTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    InitializeAODCollections(iEvent, iSetup);
    primaryVertex = vertices->ptrAt(0);
    for(auto energyScale : eventEnergyScales) {
        InitializeCandidateCollections(energyScale);
        try {
            Cutter cut(&GetAnaData().Selection("events"));
            cut(true, "events");
            ProcessEvent(cut);
        } catch(cuts::cut_failed&){}

        GetAnaData().Selection("events").fill_selection();
    }
}

void BaseTupleProducer::endJob()
{
    eventTuple.Write();
}

void BaseTupleProducer::InitializeAODCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edmEvent = &iEvent;
    triggerTools.Initialize(iEvent);
    iEvent.getByToken(electronsMiniAOD_token, pat_electrons);
    iEvent.getByToken(eleTightIdMap_token, tight_id_decisions);
    iEvent.getByToken(eleMediumIdMap_token, medium_id_decisions);
    iEvent.getByToken(eleCutBasedVetoMap_token, ele_cutBased_veto);
    iEvent.getByToken(tausMiniAOD_token, pat_taus);
    iEvent.getByToken(muonsMiniAOD_token, pat_muons);
    iEvent.getByToken(vtxMiniAOD_token, vertices);
    iEvent.getByToken(pfMETAOD_token, pfMETs);
    iEvent.getByToken(jetsMiniAOD_token, pat_jets);
    iEvent.getByToken(fatJetsMiniAOD_token, pat_fatJets);
    iEvent.getByToken(metCovMatrix_token, metCovMatrix);
    iEvent.getByToken(PUInfo_token, PUInfo);
    if(isMC) iEvent.getByToken(genWeights_token, genEvt);
    iEvent.getByToken(prunedGen_token, genParticles);
    try {iEvent.getByToken(lheEventProduct_token, lheEventProduct);} catch(...){;};

    iSetup.get<JetCorrectionsRecord>().get("AK5PF", jetCorParColl);
    jecUnc = std::shared_ptr<JetCorrectionUncertainty>(new JetCorrectionUncertainty((*jetCorParColl)["Uncertainty"]));
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
        if(tauEnergyScales.count(energyScale)) {
            const int sign = tauEnergyScales.at(energyScale);
            const double sf = 1.0 + sign * analysis::uncertainties::tau::energyUncertainty;
            const auto shiftedMomentum = tau.p4() * sf;
            tauCandidate.SetMomentum(shiftedMomentum);
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

    met = std::shared_ptr<MET>(new MET((*pfMETs)[0], *metCovMatrix));
    if(metUncertantyMap.count(energyScale)) {
        const auto shiftedMomentum = (*met)->shiftedP4(metUncertantyMap.at(energyScale));
        met->SetMomentum(shiftedMomentum);
    }
}

double BaseTupleProducer::Isolation(const pat::Electron& electron)
{
    const double sum_neutral = electron.pfIsolationVariables().sumNeutralHadronEt
                             + electron.pfIsolationVariables().sumPhotonEt
                             - 0.5 * electron.pfIsolationVariables().sumPUPt;
    const double abs_iso = electron.pfIsolationVariables().sumChargedHadronPt + std::max(sum_neutral, 0.0);
    return abs_iso / electron.pt();
}

double BaseTupleProducer::Isolation(const pat::Muon& muon)
{
    const double sum_neutral = muon.pfIsolationR04().sumNeutralHadronEt
                             + muon.pfIsolationR04().sumPhotonEt
                             - 0.5 * muon.pfIsolationR04().sumPUPt;
    const double abs_iso = muon.pfIsolationR04().sumChargedHadronPt + std::max(sum_neutral, 0.0);
    return abs_iso / muon.pt();
}

double BaseTupleProducer::Isolation(const pat::Tau& tau)
{
    return tau.tauID("byIsolationMVArun2v1DBoldDMwLTraw");
}

//  https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Recommendations_for_13_TeV_data
//  PFJetID is tuned on Uncorrected Jet values
bool BaseTupleProducer::PassPFLooseId(const pat::Jet& pat_jet)
{
    const pat::Jet& patJet = pat_jet.correctedJet("Uncorrected");
    const double abs_eta = std::abs(patJet.eta());

    if(abs_eta < 2.7 && (
        patJet.neutralHadronEnergyFraction() >= 0.99 ||
        patJet.neutralEmEnergyFraction() >= 0.99 ||
        patJet.nConstituents() <= 1)) return false;

    if(abs_eta <= 2.4 && (
        patJet.chargedHadronEnergyFraction() <= 0 ||
        patJet.chargedMultiplicity() <= 0 ||
        patJet.chargedEmEnergyFraction() >= 0.99)) return false;

    if(abs_eta > 2.7 && abs_eta <= 3.0 && (
        patJet.neutralEmEnergyFraction() >= 0.90 ||
        patJet.neutralMultiplicity() <= 2)) return false;

    if(abs_eta > 3.0 && (
        patJet.neutralEmEnergyFraction() >= 0.90 ||
        patJet.neutralMultiplicity() <= 10)) return false;

    return true;
}

void BaseTupleProducer::FillLheInfo()
{
    static constexpr int b_quark = 5;
    static const std::set<int> quarks_and_gluons = { 1, 2, 3, 4, 5, 6, 21 };
    static const std::set<int> interesting_particles = { 5, 6, 23, 24};

    if(!lheEventProduct.isValid()) {
        eventTuple().lhe_n_partons = ntuple::DefaultFillValue<UInt_t>();
        eventTuple().lhe_n_b_partons = ntuple::DefaultFillValue<UInt_t>();
        eventTuple().lhe_HT = ntuple::DefaultFillValue<Float_t>();
        return;
    }

    const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
    const std::vector<lhef::HEPEUP::FiveVector>& lheParticles = lheEvent.PUP;
    eventTuple().lhe_n_partons = 0;
    eventTuple().lhe_n_b_partons = 0;
    double HT2 = 0;
    for(size_t n = 0; n < lheParticles.size(); ++n) {
        const int absPdgId = std::abs(lheEvent.IDUP[n]);
        const int status = lheEvent.ISTUP[n];
        if(status != 1 || !quarks_and_gluons.count(absPdgId)) continue;
        eventTuple().lhe_n_partons++;
        if(absPdgId == b_quark) eventTuple().lhe_n_b_partons++;
        HT2 += std::pow(lheParticles[n][0], 2) + std::pow(lheParticles[n][1], 2);
    }
    eventTuple().lhe_HT = std::sqrt(HT2);
}

double BaseTupleProducer::GetNumberOfPileUpInteractions() const
{
    if(PUInfo.isValid()) {
        for(auto PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI) {
            if(PVI->getBunchCrossing() == 0)
                return PVI->getTrueNumInteractions();
        }
    }
    return std::numeric_limits<double>::lowest();
}

void BaseTupleProducer::ApplyBaseSelection(analysis::SelectionResultsBase& selection,
                        const std::vector<LorentzVector>& signalLeptonMomentums)
{
    selection.jets = CollectJets(signalLeptonMomentums);
    if(!runKinFit) return;

    for(size_t n = 0; n < selection.jets.size(); ++n) {
        for(size_t k = 0; k < selection.jets.size(); ++k) {
            if(k == n) continue;
            const auto& result = kinfitProducer.Fit(signalLeptonMomentums.at(0), signalLeptonMomentums.at(1),
                                                    selection.jets.at(n).GetMomentum(),
                                                    selection.jets.at(k).GetMomentum(), *met);
            selection.kinfitResults.push_back(result);
        }
    }
}

std::vector<BaseTupleProducer::ElectronCandidate> BaseTupleProducer::CollectZelectrons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectZElectron, this, _1, _2);
    return CollectObjects("Zelectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::MuonCandidate> BaseTupleProducer::CollectZmuons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectZMuon, this, _1, _2);
    return CollectObjects("Zmuons", base_selector, muons);
}

std::vector<BaseTupleProducer::ElectronCandidate> BaseTupleProducer::CollectVetoElectrons(
        const ElectronCandidate* signalElectron)
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectVetoElectron, this, _1, _2, signalElectron);
    return CollectObjects("vetoElectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::MuonCandidate> BaseTupleProducer::CollectVetoMuons(const MuonCandidate* signalMuon)
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectVetoMuon, this, _1, _2, signalMuon);
    return CollectObjects("vetoMuons", base_selector, muons);
}

std::vector<BaseTupleProducer::JetCandidate> BaseTupleProducer::CollectJets(
        const std::vector<LorentzVector>& signalLeptonMomentums)
{
    using namespace std::placeholders;
    const auto baseSelector = std::bind(&BaseTupleProducer::SelectJet, this, _1, _2, signalLeptonMomentums);
    return CollectObjects("jets", baseSelector, jets);
}

void BaseTupleProducer::SelectZElectron(const ElectronCandidate& electron, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::ZeeVeto;

    cut(true, "gt0_ele_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electron_dxy = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electron_dxy < dxy, "dxy", electron_dxy);
    const double electron_dz = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electron_dz < dz, "dz", electron_dz);
    const bool veto  = (*ele_cutBased_veto)[electron.getPtr()];
    cut(veto, "cut_based_veto");
    cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
}

void BaseTupleProducer::SelectZMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::ZmumuVeto;

    cut(true, "gt0_mu_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    cut(muon->isGlobalMuon(), "GlobalMuon");
    cut(muon->isTrackerMuon(), "trackerMuon");
    cut(muon->isPFMuon(), "PFMuon");
    cut(muon.GetIsolation() < pfRelIso04, "pfRelIso", muon.GetIsolation());
}

void BaseTupleProducer::SelectVetoElectron(const ElectronCandidate& electron, Cutter& cut,
                                           const ElectronCandidate* signalElectron) const
{
    using namespace cuts::H_tautau_2016::electronVeto;

    cut(true, "gt0_ele_cand");
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
        const auto eleMissingHits =
                electron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        cut(eleMissingHits <= missingHits, "missingHits", eleMissingHits);
        cut(electron->passConversionVeto(), "conversionVeto");
    }
    cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
    if(signalElectron) {
        const bool isNotSignal =  &(*electron) != &(*(*signalElectron));
        cut(isNotSignal, "isNotSignal");
    }
}

void BaseTupleProducer::SelectVetoMuon(const MuonCandidate& muon, Cutter& cut, const MuonCandidate* signalMuon) const
{
    using namespace cuts::H_tautau_2016::muonVeto;

    cut(true, "gt0_mu_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    cut(muon.GetIsolation() < pfRelIso04, "iso", muon.GetIsolation());

    const bool passMuonId = productionMode == ProductionMode::hh ? muon->isLooseMuon() : muon->isMediumMuon();
    cut(passMuonId, "muonID");
    if(signalMuon) {
        const bool isNotSignal =  &(*muon) != &(*(*signalMuon));
        cut(isNotSignal, "isNotSignal");
    }
}

void BaseTupleProducer::SelectJet(const JetCandidate& jet, Cutter& cut,
                                  const std::vector<LorentzVector>& signalLeptonMomentums) const
{
    using namespace cuts::H_tautau_2016::jetID;

    cut(true, "gt0_jet_cand");
    const LorentzVector& p4 = jet.GetMomentum();
    cut(p4.Pt() > pt, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    cut(PassPFLooseId(*jet), "jet_id");
    for(size_t n = 0; n < signalLeptonMomentums.size(); ++n) {
        std::ostringstream cut_name;
        cut_name << "deltaR_lep" << n + 1;
        const double deltaR = ROOT::Math::VectorUtil::DeltaR(p4, signalLeptonMomentums.at(n));
        cut(deltaR > deltaR_signalObjects, cut_name.str(), deltaR);
    }
}

void BaseTupleProducer::FillEventTuple(const analysis::SelectionResultsBase& selection)
{
    using namespace analysis;
    static constexpr float default_value = ntuple::DefaultFillValue<Float_t>();
    static constexpr int default_int_value = ntuple::DefaultFillValue<Int_t>();

    eventTuple().run  = edmEvent->id().run();
    eventTuple().lumi = edmEvent->id().luminosityBlock();
    eventTuple().evt  = edmEvent->id().event();
    eventTuple().eventEnergyScale = static_cast<int>(eventEnergyScale);

    eventTuple().weightevt = isMC ? genEvt->weight() : 1;
    FillLheInfo();

    eventTuple().npv = vertices->size();
    eventTuple().npu = GetNumberOfPileUpInteractions();

    // HTT candidate
    eventTuple().SVfit_p4 = selection.svfitResult.momentum;

    // Leg 2, tau
    const TauCandidate& tau = selection.GetSecondLeg();
    eventTuple().p4_2     = analysis::LorentzVectorM(tau.GetMomentum());
    eventTuple().q_2      = tau.GetCharge();
    eventTuple().d0_2     = Calculate_dxy(tau->vertex(), primaryVertex->position(), tau.GetMomentum());
    eventTuple().dZ_2     = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get())->dz();
    eventTuple().iso_2    = tau.GetIsolation();
    eventTuple().id_e_mva_nt_loose_1 = default_value;
    eventTuple().gen_match_2 = isMC ? gen_truth::genMatch(tau->p4(), *genParticles) : default_int_value;

    eventTuple().tauIDs_2.insert(tau->tauIDs().begin(), tau->tauIDs().end());

    // MET
    eventTuple().pfMET_p4 = met->GetMomentum();
    eventTuple().pfMET_cov = met->GetCovMatrix();

    for(const JetCandidate& jet : selection.jets) {
        const LorentzVector& p4 = jet.GetMomentum();
        eventTuple().jets_p4.push_back(p4);
        eventTuple().jets_rawf.push_back((jet->correctedJet("Uncorrected").pt() ) / p4.Pt());
        eventTuple().jets_mva.push_back(jet->userFloat("pileupJetId:fullDiscriminant"));
        eventTuple().jets_csv.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        //eventTuple().jets_partonFlavour.push_back(jet->partonFlavour());
        eventTuple().jets_hadronFlavour.push_back(jet->hadronFlavour());
    }

    for(const JetCandidate& jet : fatJets) {
        const LorentzVector& p4 = jet.GetMomentum();
        eventTuple().fatJets_p4.push_back(p4);
        eventTuple().fatJets_csv.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        eventTuple().fatJets_m_pruned.push_back(GetUserFloat(jet, "ak8PFJetsCHSPrunedMass"));
        eventTuple().fatJets_m_filtered.push_back(GetUserFloat(jet, "ak8PFJetsCHSFilteredMass"));
        eventTuple().fatJets_m_trimmed.push_back(GetUserFloat(jet, "ak8PFJetsCHSTrimmedMass"));
        eventTuple().fatJets_m_softDrop.push_back(GetUserFloat(jet, "ak8PFJetsCHSSoftDropMass"));
        eventTuple().fatJets_n_subjettiness_tau1.push_back(GetUserFloat(jet, "NjettinessAK8:tau1"));
        eventTuple().fatJets_n_subjettiness_tau2.push_back(GetUserFloat(jet, "NjettinessAK8:tau2"));
        eventTuple().fatJets_n_subjettiness_tau3.push_back(GetUserFloat(jet, "NjettinessAK8:tau3"));

        if(!jet->hasSubjets("SoftDrop")) continue;
        const size_t parentIndex = eventTuple().fatJets_p4.size() - 1;
        const auto& sub_jets = jet->subjets("SoftDrop");
        for(const auto& sub_jet : sub_jets) {
            eventTuple().subJets_p4.push_back(analysis::LorentzVector(sub_jet->p4()));
            eventTuple().subJets_csv.push_back(sub_jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            eventTuple().subJets_parentIndex.push_back(parentIndex);
        }
    }

    for(const kin_fit::FitResults& result : selection.kinfitResults) {
        eventTuple().kinFit_m.push_back(result.mass);
        eventTuple().kinFit_chi2.push_back(result.chi2);
        eventTuple().kinFit_convergence.push_back(result.convergence);
    }

    eventTuple().dilepton_veto  = selection.Zveto;
    eventTuple().extraelec_veto = selection.electronVeto;
    eventTuple().extramuon_veto = selection.muonVeto;
}
