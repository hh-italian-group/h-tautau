#include "h-tautau/Production/interface/BaseTupleProducer.h"
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
    metCovMatrix_token(consumes<MetCovMatrix>(iConfig.getParameter<edm::InputTag>("metCov"))),
    PUInfo_token(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
    lheEventProduct_token(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
    genWeights_token(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
    prunedGen_token(consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
    isMC(iConfig.getParameter<bool>("isMC")),
    applyTriggerMatch(iConfig.getParameter<bool>("applyTriggerMatch")),
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
}

void BaseTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using analysis::EventEnergyScale;
    static const std::set<EventEnergyScale> energyScales = {
        EventEnergyScale::Central, EventEnergyScale::TauUp, EventEnergyScale::TauDown,
        EventEnergyScale::JetUp, EventEnergyScale::JetDown
    };

    InitializeAODCollections(iEvent, iSetup);
    primaryVertex = vertices->ptrAt(0);
    for(auto energyScale : energyScales) {
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
            const double sf = 1.0 + sign * cuts::Htautau_2015::tauCorrections::energyUncertainty;
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
    const double sum_neutral = muon.pfIsolationR03().sumNeutralHadronEt
                             + muon.pfIsolationR03().sumPhotonEt
                             - 0.5 * muon.pfIsolationR03().sumPUPt;
    const double abs_iso = muon.pfIsolationR03().sumChargedHadronPt + std::max(sum_neutral, 0.0);
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
    //TLorentzVector momentum = jet->GetMomentum();
    const pat::Jet& patJet = pat_jet.correctedJet("Uncorrected");
    //momentum.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), jet.mass());
    if(std::abs(patJet.eta())<3.0) {
        //if(momentum.E() == 0)                                  return false;
        if(patJet.neutralHadronEnergyFraction() > 0.99)   return false;
        if(patJet.neutralEmEnergyFraction()     > 0.99)   return false;
        if(patJet.nConstituents() <  1)                   return false;
        if(patJet.chargedHadronEnergyFraction() <= 0 && std::abs(patJet.eta()) < 2.4 ) return false;
        if(patJet.chargedEmEnergyFraction() >  0.99  && std::abs(patJet.eta()) < 2.4 ) return false;
        if(patJet.chargedMultiplicity()     <= 0      && std::abs(patJet.eta()) < 2.4 ) return false;
    }
    if(std::abs(patJet.eta())>3.0) {
        if(patJet.neutralEmEnergyFraction()     > 0.90)   return false;
        if(patJet.neutralMultiplicity() < 10 )            return false;
    }
    return true;
}

std::pair<double,int> BaseTupleProducer::ComputeHtValue()
{
    if(!lheEventProduct.isValid())
        return std::make_pair(-999.99,-1);

    const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
    double lheHt = 0.;
    int lheNOutPartons = 0;
    size_t numParticles = lheParticles.size();
    for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
        int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
        int status = lheEvent.ISTUP[idxParticle];
        if ( status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
            lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.)
                    + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
            ++lheNOutPartons;
        }
    }
    return std::make_pair(lheHt,lheNOutPartons);
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
    selection.bjets = CollectBJets(signalLeptonMomentums);

    for(size_t n = 0; n < selection.jets.size(); ++n) {
        for(size_t k = n + 1; k < selection.jets.size(); ++k) {
            const std::vector<LorentzVector> jet_momentums = {
                selection.jets.at(n).GetMomentum(), selection.jets.at(k).GetMomentum()
            };
            const auto& result = kinfitProducer.Fit(signalLeptonMomentums, jet_momentums, *met);
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

std::vector<BaseTupleProducer::JetCandidate> BaseTupleProducer::CollectBJets(
        const std::vector<LorentzVector>& signalLeptonMomentums)
{
    using namespace std::placeholders;
    static const std::string& btagDiscName = "pfCombinedInclusiveSecondaryVertexV2BJetTags";
    const auto jetCsvOrdering = [&](const JetCandidate& jet1, const JetCandidate& jet2) -> bool {
        return jet1->bDiscriminator(btagDiscName) > jet2->bDiscriminator(btagDiscName);
    };
    const auto baseSelector = std::bind(&BaseTupleProducer::SelectBJet, this, _1, _2, signalLeptonMomentums);
    return CollectObjects("bjets", baseSelector, jets, jetCsvOrdering);
}

void BaseTupleProducer::SelectZElectron(const ElectronCandidate& electron, Cutter& cut) const
{
    using namespace cuts::Htautau_2015::ETau;

    cut(true, "gt0_ele_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    cut(p4.pt() > ZeeVeto::pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < ZeeVeto::eta, "eta", p4.eta());
    const double electronD0 = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electronD0 < ZeeVeto::d0, "dxy", electronD0);
    const double electronDZ = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electronDZ < ZeeVeto::dz, "dz", electronDZ);
    const bool veto  = (*ele_cutBased_veto)[electron.getPtr()];
    cut(veto, "cut_based_veto");
    cut(electron.GetIsolation() < ZeeVeto::pfRelIso, "iso", electron.GetIsolation());
}

void BaseTupleProducer::SelectZMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::Htautau_2015::MuTau;

    cut(true, "gt0_mu_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > ZmumuVeto::pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < ZmumuVeto::eta, "eta", p4.eta());
    const double muonDZ = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muonDZ < muonID::dz, "dz", muonDZ);
    const double muonDB = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muonDB < muonID::dB, "dxy", muonDB);
    cut(muon->isTrackerMuon(), "trackerMuon");
    cut(muon->isGlobalMuon(), "GlobalMuon");
    cut(muon->isPFMuon(), "PFMuon");
    cut(muon.GetIsolation() < 0.3, "pFRelIso", muon.GetIsolation());
}

void BaseTupleProducer::SelectVetoElectron(const ElectronCandidate& electron, Cutter& cut,
                                           const ElectronCandidate* signalElectron) const
{
    using namespace cuts::Htautau_2015;

    cut(true, "gt0_ele_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    cut(p4.pt() > electronVeto::pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < electronVeto::eta_high, "eta", p4.eta());
    const double electronD0 = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electronD0 < electronVeto::d0, "dxy", electronD0);
    const double electronDZ = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electronDZ < electronVeto::dz, "dz", electronDZ);
    const bool isMedium  = (*medium_id_decisions)[electron.getPtr()];
    cut(isMedium, "electronMVAMediumID");
    const int eleMissingHits = electron->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
    cut(eleMissingHits <= electronVeto::missingHits, "missingHits", eleMissingHits);
    cut(electron->passConversionVeto(),"conversionVeto");
    cut(electron.GetIsolation() < electronVeto::pFRelIso, "iso", electron.GetIsolation());
    if(signalElectron) {
        const double deltaR = ROOT::Math::VectorUtil::DeltaR(p4, signalElectron->GetMomentum());
        cut(deltaR > 0.05, "deltaR_signal", deltaR);
    }
}

void BaseTupleProducer::SelectVetoMuon(const MuonCandidate& muon, Cutter& cut, const MuonCandidate* signalMuon) const
{
    using namespace cuts::Htautau_2015;

    cut(true, "gt0_mu_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > muonVeto::pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < muonVeto::eta, "eta", p4.eta());
    const double muonDB = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muonDB < muonVeto::dB, "dxy", muonDB);
    const double muonDZ = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muonDZ < muonVeto::dz, "dz", muonDZ);
    cut(muon.GetIsolation() < muonVeto::pfRelIso, "iso", muon.GetIsolation());
    cut(muon->isMediumMuon(), "muonID");
    if(signalMuon) {
        const double deltaR = ROOT::Math::VectorUtil::DeltaR(p4, signalMuon->GetMomentum());
        cut(deltaR > 0.05, "deltaR_signal", deltaR);
    }
}

void BaseTupleProducer::SelectJet(const JetCandidate& jet, Cutter& cut,
                                  const std::vector<LorentzVector>& signalLeptonMomentums) const
{
    using namespace cuts::Htautau_2015;
    using namespace cuts::Htautau_2015::jetID;

    cut(true, "gt0_jet_cand");
    const LorentzVector& p4 = jet.GetMomentum();
    cut(p4.Pt() > pt_loose, "pt_loose", p4.Pt());
    cut(std::abs( p4.Eta() ) < eta, "eta", p4.Eta());
    const bool jetPFID = PassPFLooseId(*jet);
    cut(jetPFID, "jet_id");
    for(size_t n = 0; n < signalLeptonMomentums.size(); ++n) {
        std::ostringstream cut_name;
        cut_name << "deltaR_lep" << n + 1;
        const double deltaR = ROOT::Math::VectorUtil::DeltaR(p4, signalLeptonMomentums.at(n));
        cut(deltaR > deltaR_signalObjects, cut_name.str(), deltaR);
    }
}

void BaseTupleProducer::SelectBJet(const JetCandidate& jet, Cutter& cut,
                                   const std::vector<LorentzVector>& signalLeptonMomentums) const
{
    using namespace cuts::Htautau_2015;
    using namespace cuts::Htautau_2015::jetID;

    SelectJet(jet, cut, signalLeptonMomentums);
    const LorentzVector& p4 = jet.GetMomentum();
    cut(std::abs(p4.eta()) < btag::eta, "eta", p4.eta());
    const double csvValue = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    cut(csvValue > btag::CSV, "btaggin", csvValue);
}

void BaseTupleProducer::FillEventTuple(const analysis::SelectionResultsBase& selection)
{
    using namespace analysis;
    static const float default_value = ntuple::DefaultFillValue<Float_t>();
    static const std::set<int> selected_lhe_pdgs = { 1, 2, 3, 4, 5, 6, 11, 13, 15, 21 };

    eventTuple().run  = edmEvent->id().run();
    eventTuple().lumi = edmEvent->id().luminosityBlock();
    eventTuple().evt  = edmEvent->id().event();
    eventTuple().eventEnergyScale = static_cast<int>(eventEnergyScale);

    eventTuple().weightevt = isMC ? genEvt->weight() : default_value;
    if(lheEventProduct.isValid()) {
        const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
        const std::vector<lhef::HEPEUP::FiveVector>& lheParticles = lheEvent.PUP;
        for(size_t n = 0; n < lheParticles.size(); ++n) {
            const int pdg_id = lheEvent.IDUP.at(n);
            const int status = lheEvent.ISTUP.at(n);
            if(status != 1 || !selected_lhe_pdgs.count(std::abs(pdg_id))) continue;
            const auto& momentum = lheParticles.at(n);
            const analysis::LorentzVectorXYZ p4(momentum[0], momentum[1], momentum[2], momentum[3]);
            eventTuple().lhe_particle_pdg.push_back(pdg_id);
            eventTuple().lhe_particle_p4.push_back(analysis::LorentzVectorM(p4));
        }
    }

    eventTuple().npv = vertices->size();
    eventTuple().npu = GetNumberOfPileUpInteractions();

    // HTT candidate
    eventTuple().SVfit_p4 = selection.svfitResult.momentum;

    // Leg 2, tau
    const TauCandidate& tau = selection.GetSecondLeg();
    eventTuple().p4_2     = analysis::LorentzVectorM(tau.GetMomentum());
    eventTuple().q_2      = tau.GetCharge();
    eventTuple().pfmt_2   = Calculate_MT(tau.GetMomentum(), met->GetMomentum().Pt(), met->GetMomentum().Phi());
    eventTuple().d0_2     = Calculate_dxy(tau->vertex(), primaryVertex->position(), tau.GetMomentum());
    eventTuple().dZ_2     = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get())->dz();
    eventTuple().iso_2    = tau.GetIsolation();
    eventTuple().id_e_mva_nt_loose_1 = default_value;
    eventTuple().gen_match_2 = isMC ? gen_truth::genMatch(tau->p4(), *genParticles) : default_value;

    eventTuple().againstElectronLooseMVA6_2   = tau->tauID("againstElectronLooseMVA6");
    eventTuple().againstElectronMediumMVA6_2  = tau->tauID("againstElectronMediumMVA6");
    eventTuple().againstElectronTightMVA6_2   = tau->tauID("againstElectronTightMVA6");
    eventTuple().againstElectronVLooseMVA6_2  = tau->tauID("againstElectronVLooseMVA6");
    eventTuple().againstElectronVTightMVA6_2  = tau->tauID("againstElectronVTightMVA6");

    eventTuple().againstMuonLoose3_2          = tau->tauID("againstMuonLoose3");
    eventTuple().againstMuonTight3_2          = tau->tauID("againstMuonTight3");

    eventTuple().byCombinedIsolationDeltaBetaCorrRaw3Hits_2 = tau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
    eventTuple().byIsolationMVA3newDMwLTraw_2               = tau->tauID("byIsolationMVArun2v1DBnewDMwLTraw");
    eventTuple().byIsolationMVA3oldDMwLTraw_2               = tau->tauID("byIsolationMVArun2v1DBoldDMwLTraw");
    eventTuple().byIsolationMVA3newDMwoLTraw_2              = default_value;
    eventTuple().byIsolationMVA3oldDMwoLTraw_2              = default_value;

    eventTuple().byVLooseIsolationMVArun2v1DBoldDMwLT_2     = tau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byLooseIsolationMVArun2v1DBoldDMwLT_2      = tau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byMediumIsolationMVArun2v1DBoldDMwLT_2     = tau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byTightIsolationMVArun2v1DBoldDMwLT_2      = tau->tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    eventTuple().byVTightIsolationMVArun2v1DBoldDMwLT_2     = tau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

    eventTuple().decayModeFindingOldDMs_2 = tau->tauID("decayModeFinding");

    // MET
    eventTuple().pfMET_p4 = met->GetMomentum();
    eventTuple().pfMET_cov = met->GetCovMatrix();

    for(const JetCandidate& jet : selection.jets){
        const LorentzVector& p4 = jet.GetMomentum();
        eventTuple().jets_p4.push_back(p4);
        eventTuple().jets_rawf.push_back((jet->correctedJet("Uncorrected").pt() ) / p4.Pt());
        eventTuple().jets_mva.push_back(jet->userFloat("pileupJetId:fullDiscriminant"));
        eventTuple().jets_csv.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
        eventTuple().jets_partonFlavour .push_back(jet->partonFlavour());
    }

    for(const kin_fit::FitResults& result : selection.kinfitResults) {
        eventTuple().kinFit_m.push_back(result.mass);
        eventTuple().kinFit_chi2.push_back(result.chi2);
        eventTuple().kinFit_probability.push_back(result.probability);
        eventTuple().kinFit_convergence.push_back(result.convergence);
    }

    eventTuple().dilepton_veto  = selection.Zveto;
    eventTuple().extraelec_veto = selection.electronVeto;
    eventTuple().extramuon_veto = selection.muonVeto;
}
