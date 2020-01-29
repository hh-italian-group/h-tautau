#include "TROOT.h"
#include "h-tautau/Production/interface/BaseTupleProducer.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "../interface/GenTruthTools.h"
#include "h-tautau/Analysis/include/MetFilters.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/Analysis/include/EventInfo.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "h-tautau/Core/include/DiscriminatorIdResults.h"


int TupleStore::tuple_counter = 0;

std::shared_ptr<ntuple::EventTuple> TupleStore::eventTuple_ptr;

ntuple::EventTuple& TupleStore::GetTuple()
{
  if(tuple_counter == 0){
    TFile& file = edm::Service<TFileService>()->file();
    file.SetCompressionAlgorithm(ROOT::kLZ4);
    file.SetCompressionLevel(4);
    eventTuple_ptr = ntuple::CreateEventTuple("events",&file,false,ntuple::TreeState::Full);
  }
  ++tuple_counter;
  return *eventTuple_ptr;
}

void TupleStore::ReleaseEventTuple()
{
  if(tuple_counter==0)
    throw analysis::exception("Tuple Counter equal zero.");
  --tuple_counter;
  if(tuple_counter == 0){
    eventTuple_ptr->Write();
    eventTuple_ptr.reset();
  }
}

BaseTupleProducer::BaseTupleProducer(const edm::ParameterSet& iConfig, analysis::Channel _channel) :
    treeName(ToString(_channel)),
    anaData(&edm::Service<TFileService>()->file(), treeName + "_stat"),
    electronsMiniAOD_token(mayConsume<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronSrc"))),
    tausMiniAOD_token(mayConsume<std::vector<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
    muonsMiniAOD_token(mayConsume<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    vtxMiniAOD_token(mayConsume<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
    pfMETAOD_token(mayConsume<edm::View<pat::MET> >(iConfig.getParameter<edm::InputTag>("pfMETSrc"))),
    genMETAOD_token(mayConsume<edm::View<reco::GenMET> >(iConfig.getParameter<edm::InputTag>("genMetSrc"))),
    jetsMiniAOD_token(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jetSrc"))),
    fatJetsMiniAOD_token(mayConsume<std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("fatJetSrc"))),
    PUInfo_token(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PUInfo"))),
    lheEventProduct_token(mayConsume<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEventProducts"))),
    genWeights_token(mayConsume<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEventInfoProduct"))),
    topGenEvent_token(mayConsume<TtGenEvent>(iConfig.getParameter<edm::InputTag>("topGenEvent"))),
    genParticles_token(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    genJets_token(mayConsume<edm::View<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    period(analysis::EnumNameMap<analysis::Period>::GetDefault().Parse(
                       iConfig.getParameter<std::string>("period"))),
    isMC(iConfig.getParameter<bool>("isMC")),
    applyTriggerMatch(iConfig.getParameter<bool>("applyTriggerMatch")),
    applyTriggerMatchCut(iConfig.getParameter<bool>("applyTriggerMatchCut")),
    runSVfit(iConfig.getParameter<bool>("runSVfit")),
    applyTriggerCut(iConfig.getParameter<bool>("applyTriggerCut")),
    storeLHEinfo(iConfig.getParameter<bool>("storeLHEinfo")),
    applyRecoilCorr(iConfig.getParameter<bool>("applyRecoilCorr")),
    nJetsRecoilCorr(iConfig.getParameter<int>("nJetsRecoilCorr")),
    saveGenTopInfo(iConfig.getParameter<bool>("saveGenTopInfo")),
    saveGenBosonInfo(iConfig.getParameter<bool>("saveGenBosonInfo")),
    saveGenJetInfo(iConfig.getParameter<bool>("saveGenJetInfo")),
    saveGenParticleInfo(iConfig.getParameter<bool>("saveGenParticleInfo")),
    isEmbedded(iConfig.getParameter<bool>("isEmbedded")),
    //eventTuple_ptr(ntuple::CreateEventTuple(ToString(_channel),&edm::Service<TFileService>()->file(),false,ntuple::TreeState::Full)),
    eventTuple(TupleStore::GetTuple()),
    //eventTuple(*eventTuple_ptr),
    triggerTools(mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "SIM")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "PAT")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "SIMembedding")),
                 mayConsume<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "MERGE")),
                 consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales")),
                 consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects")),
                 mayConsume<BXVector<l1t::Tau>>(edm::InputTag("caloStage2Digis", "Tau", isEmbedded ? "SIMembedding" : "RECO")),
                 edm::FileInPath(iConfig.getParameter<std::string>("triggerCfg")).fullPath(),
                 _channel, isEmbedded)
{
    root_ext::HistogramFactory<TH1D>::LoadConfig(
            edm::FileInPath("h-tautau/Production/data/histograms.cfg").fullPath());

    // if(period == analysis::Period::Run2016){
    //     badPFMuonFilter_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("badPFMuonFilter"));
    //     badChCandidateFilter_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("badChCandidateFilter"));
    // }

    m_rho_token = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));

    if(runSVfit)
//        svfitProducer = std::shared_ptr<analysis::sv_fit::FitProducer>(new analysis::sv_fit::FitProducer(
//            edm::FileInPath("TauAnalysis/SVfitStandalone/data/svFitVisMassAndPtResolutionPDF.root").fullPath()));
        svfitProducer = std::shared_ptr<analysis::sv_fit::FitProducer>(new analysis::sv_fit::FitProducer());

    if(applyRecoilCorr)
        recoilPFMetCorrector = std::shared_ptr<RecoilCorrector>(new RecoilCorrector(
            edm::FileInPath("HTT-utilities/RecoilCorrections/data/TypeIPFMET_2016BCD.root").fullPath()));

    const edm::ParameterSet& customMetFilters = iConfig.getParameterSet("customMetFilters");
    for(const auto& filterName : customMetFilters.getParameterNames()) {
            customMetFilters_token[filterName] =
                mayConsume<bool>(customMetFilters.getParameter<edm::InputTag>(filterName));
        }
}

void BaseTupleProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    InitializeAODCollections(iEvent, iSetup);
    primaryVertex = vertices->ptrAt(0);
    InitializeCandidateCollections();
    try {
        Cutter cut(&GetAnaData().Selection());
        cut(true,"events");
        ProcessEvent(cut);
    } catch(cuts::cut_failed&){}

    GetAnaData().Selection().fill_selection();

}

void BaseTupleProducer::endJob()
{
    //eventTuple.Write();
    TupleStore::ReleaseEventTuple();
}

void BaseTupleProducer::InitializeAODCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edmEvent = &iEvent;
    eventId = iEvent.id();
    triggerTools.Initialize(iEvent,!isMC && !isEmbedded);

    iEvent.getByToken(electronsMiniAOD_token, pat_electrons);
    iEvent.getByToken(tausMiniAOD_token, pat_taus);
    iEvent.getByToken(muonsMiniAOD_token, pat_muons);
    iEvent.getByToken(vtxMiniAOD_token, vertices);
    iEvent.getByToken(pfMETAOD_token, pfMETs);
    iEvent.getByToken(jetsMiniAOD_token, pat_jets);
    iEvent.getByToken(fatJetsMiniAOD_token, pat_fatJets);
    iEvent.getByToken(PUInfo_token, PUInfo);
    if(isMC) {
        if(!isEmbedded)
          iEvent.getByToken(genMETAOD_token, genMET);
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

void BaseTupleProducer::InitializeCandidateCollections()
{
    using METUncertainty = pat::MET::METUncertainty;
    electrons.clear();
    for(size_t n = 0; n < pat_electrons->size(); ++n) {
        const edm::Ptr<pat::Electron> ele_ptr(pat_electrons, n);
        electrons.push_back(ElectronCandidate(ele_ptr, Isolation(*ele_ptr)));
    }

    muons.clear();
    for(const auto& muon : *pat_muons)
        muons.push_back(MuonCandidate(muon, Isolation(muon)));

    met = std::shared_ptr<MET>(new MET((*pfMETs)[0], (*pfMETs)[0].getSignificanceMatrix()));

    taus.clear();
    for(const auto& tau : *pat_taus) {
        TauCandidate tauCandidate(tau,0);
        taus.push_back(tauCandidate);
    }

    jets.clear();
    for(const auto& jet : *pat_jets) {
        JetCandidate jetCandidate(jet);
        jets.push_back(jetCandidate);
    }

    fatJets.clear();
    for(const auto& jet : * pat_fatJets)
        fatJets.push_back(JetCandidate(jet));
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

    if(period == analysis::Period::Run2018)
    {
        if(abs_eta <= 2.6 && (
                         patJet.neutralHadronEnergyFraction() >= 0.9 ||
                         patJet.neutralEmEnergyFraction() >= 0.9 ||
                         patJet.nConstituents() <= 1 ||
                         patJet.chargedHadronEnergyFraction() <= 0 ||
                         patJet.chargedMultiplicity() <= 0)) return false;
        if(abs_eta > 2.6 && abs_eta <= 2.7 && (
                         patJet.neutralHadronEnergyFraction() >= 0.9 ||
                         patJet.neutralEmEnergyFraction() >= 0.99 ||
                         patJet.chargedMultiplicity() <= 0)) return false;

        if(abs_eta > 2.7 && abs_eta <= 3.0 && (
                                           patJet.neutralEmEnergyFraction() <= 0.02 ||
                                           patJet.neutralEmEnergyFraction() >= 0.99 ||
                                           patJet.neutralMultiplicity() <= 2)) return false;

        if(abs_eta > 3.0 && (
                         patJet.neutralEmEnergyFraction() >= 0.9 ||
                         patJet.neutralHadronEnergyFraction() <= 0.2 ||
                         patJet.neutralMultiplicity() <= 10)) return false;

    }
    return true;
}

void BaseTupleProducer::FillLheInfo()
{
    if(!lheEventProduct.isValid()) {
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
    if(storeLHEinfo){
    	for(size_t n = 0; n < lheSummary.index.size(); ++n){
         	eventTuple().lhe_index.push_back(lheSummary.index.at(n));
        	eventTuple().lhe_pdgId.push_back(lheSummary.pdgId.at(n));
        	eventTuple().lhe_first_mother_index.push_back(lheSummary.first_mother_index.at(n));
            eventTuple().lhe_last_mother_index.push_back(lheSummary.last_mother_index.at(n));
        	eventTuple().lhe_p4.push_back(ntuple::LorentzVectorM(lheSummary.p4.at(n)));
    	}
    }

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
    static constexpr int electronNeutrinoPdgId = 12, muonNeutrinoPdgId = 14, tauNeutrinoPdgId = 16;
    static const std::set<int> chargedLeptons = { electronPdgId, muonPdgId, tauPdgId };
    static const std::set<int> neutralLeptons = { electronNeutrinoPdgId, muonNeutrinoPdgId, tauNeutrinoPdgId };
    static const std::set<int> bosons = { 23, 24, 25, 35 };

    std::vector<const reco::GenParticle*> particles_to_store;
    std::vector<const reco::GenParticle*> particles_to_store_backup;

    std::map<int, size_t> particle_counts;
    for(const auto& particle : *genParticles) {
        const auto& flag = particle.statusFlags();
        if(!flag.isPrompt() || !flag.isLastCopy() || !flag.fromHardProcess()) continue;
        const int abs_pdg = std::abs(particle.pdgId());
        ++particle_counts[abs_pdg];
        if(saveGenBosonInfo) {
            if(bosons.count(abs_pdg))
                particles_to_store.push_back(&particle);
            else if(chargedLeptons.count(abs_pdg) || neutralLeptons.count(abs_pdg))
                particles_to_store_backup.push_back(&particle);
        }
    }

    if(saveGenBosonInfo && particles_to_store.empty()){
        if(particles_to_store_backup.size() < 2)
            throw analysis::exception("Less than 2 prompt leptons in the event");
        particles_to_store = particles_to_store_backup;
    }

    if(saveGenBosonInfo && particles_to_store.empty()) {
        throw analysis::exception("Particles to store is empty for event %1%.")
                % analysis::EventIdentifier(edmEvent->id().run(), edmEvent->id().luminosityBlock(),
                                            edmEvent->id().event());
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

    auto returnIndex = [&](const reco::GenParticle* particle_ptr){

		int particle_index = -1;
		if(particle_ptr != nullptr){
            particle_index = static_cast<int>(particle_ptr - genParticles->data());
		    if(particle_index > static_cast<int>(genParticles->size()) || particle_index < 0) {
                if(saveGenTopInfo && std::abs(particle_ptr->pdgId()) == 6) {
                    particle_index = -10;
                } else {
		                  throw analysis::exception("Particle index = %1% for particle with pdgId = %2% exceeds"
                                                    " the size of gen particles collection = %3%.")
                                                    % particle_index % particle_ptr->pdgId() % genParticles->size();
                }
            }
		}
		return particle_index;
	};

    auto fillGenInfo = [&](const reco::GenParticle* particle) {
    	int index = returnIndex(particle);
        eventTuple().genParticles_index.push_back(index);
    	eventTuple().genParticles_vertex.push_back(ntuple::Point3D(particle->vertex()));
    	eventTuple().genParticles_pdg.push_back(particle->pdgId());
    	eventTuple().genParticles_status.push_back(particle->status());
    	eventTuple().genParticles_statusFlags.push_back(static_cast<uint16_t>(
            particle->statusFlags().flags_.to_ulong()));
        eventTuple().genParticles_p4.push_back(ntuple::LorentzVectorM(particle->p4()));

        if(index >= 0) {
        	for(size_t mother_id = 0; mother_id < particle->numberOfMothers(); ++mother_id) {
                eventTuple().genParticles_rel_pIndex.push_back(index);
        		const auto mother_ptr = dynamic_cast<const reco::GenParticle*>(particle->mother(mother_id));
        		int mother_index = returnIndex(mother_ptr);
        		eventTuple().genParticles_rel_mIndex.push_back(mother_index);
        	}
        }

    };



    if(saveGenParticleInfo){
    	for(const auto& particle : *genParticles) {
    	    fillGenInfo(&particle);
    	}
	}
    else{
    	for(const auto& particle : particles_to_store) {
    	    fillGenInfo(particle);
    	}
    }


}

void BaseTupleProducer::FillGenJetInfo()
{
    static constexpr int b_flavour = 5, c_flavour = 4;
    static constexpr double pt_cut = 5;
    eventTuple().genJets_nTotal = genJets->size();

    std::map<int, size_t> hf_counts;

    for(const JetCandidate& jet : jets) {
        ++hf_counts[std::abs(jet->hadronFlavour())];
    }
    eventTuple().jets_nTotal_hadronFlavour_b = hf_counts[b_flavour];
    eventTuple().jets_nTotal_hadronFlavour_c = hf_counts[c_flavour];

    if(!saveGenJetInfo) return;

    for(const reco::GenJet& gen_jet : *genJets) {
        if(gen_jet.pt() <= pt_cut) continue;
        eventTuple().genJets_p4.push_back(ntuple::LorentzVectorE(gen_jet.p4()));

        const auto findRecoJetFlavour = [&]() {
            for(const JetCandidate& reco_jet : jets) {
                if(reco_jet->genJet() == &gen_jet)
                    return reco_jet->hadronFlavour();
            }
            return ntuple::DefaultFillValue<int>();
        };

        const auto flavour = findRecoJetFlavour();
        eventTuple().genJets_hadronFlavour.push_back(flavour);
    }
}

void BaseTupleProducer::FillOtherLeptons(const std::vector<ElectronCandidate>& other_electrons,
                                         const std::vector<MuonCandidate>& other_muons)
{
    for (const auto electron : other_electrons){
        eventTuple().other_lepton_p4.push_back(ntuple::LorentzVectorM(electron.GetMomentum()));
        eventTuple().other_lepton_q.push_back(electron.GetCharge());
        eventTuple().other_lepton_type.push_back(static_cast<int>(analysis::LegType::e));
        eventTuple().other_lepton_iso.push_back(electron.GetIsolation());
        eventTuple().other_lepton_elePassConversionVeto.push_back(electron->passConversionVeto());
        analysis::DiscriminatorIdResults eleId_iso;
        eleId_iso.SetResult(analysis::DiscriminatorWP::Loose,electron->electronID("mvaEleID-Fall17-iso-V1-wpLoose") > 0.5);
        eleId_iso.SetResult(analysis::DiscriminatorWP::Medium,electron->electronID("mvaEleID-Fall17-iso-V1-wp90") > 0.5);
        eleId_iso.SetResult(analysis::DiscriminatorWP::Tight,electron->electronID("mvaEleID-Fall17-iso-V1-wp80") > 0.5);
        eventTuple().other_lepton_eleId_iso.push_back(eleId_iso.GetResultBits());
        analysis::DiscriminatorIdResults eleId_noIso;
        eleId_noIso.SetResult(analysis::DiscriminatorWP::Loose,electron->electronID("mvaEleID-Fall17-noIso-V1-wpLoose") > 0.5);
        eleId_noIso.SetResult(analysis::DiscriminatorWP::Medium,electron->electronID("mvaEleID-Fall17-noIso-V1-wp90") > 0.5);
        eleId_noIso.SetResult(analysis::DiscriminatorWP::Tight,electron->electronID("mvaEleID-Fall17-noIso-V1-wp80") > 0.5);
        eventTuple().other_lepton_eleId_noIso.push_back(eleId_noIso.GetResultBits());
        eventTuple().other_lepton_muonId.push_back(0);
        if(isMC) {
            const auto match = analysis::gen_truth::LeptonGenMatch(analysis::LorentzVectorM(electron.GetMomentum()), *genParticles);
            eventTuple().other_lepton_gen_match.push_back(static_cast<int>(match.match));
            const auto matched_p4 = match.gen_particle ? match.gen_particle->p4() : analysis::LorentzVectorXYZ();
            eventTuple().other_lepton_gen_p4.push_back(ntuple::LorentzVectorM(matched_p4));
        }
    }

    for (const auto muon : other_muons){
        eventTuple().other_lepton_p4.push_back(ntuple::LorentzVectorM(muon.GetMomentum()));
        eventTuple().other_lepton_q.push_back(muon.GetCharge());
        eventTuple().other_lepton_type.push_back(static_cast<int>(analysis::LegType::mu));
        eventTuple().other_lepton_iso.push_back(muon.GetIsolation());
        eventTuple().other_lepton_elePassConversionVeto.push_back(false);
        eventTuple().other_lepton_eleId_iso.push_back(0);
        eventTuple().other_lepton_eleId_noIso.push_back(0);
        analysis::DiscriminatorIdResults muonId;
        muonId.SetResult(analysis::DiscriminatorWP::Loose,muon->isLooseMuon());
        muonId.SetResult(analysis::DiscriminatorWP::Medium,muon->isMediumMuon());
        muonId.SetResult(analysis::DiscriminatorWP::Tight,muon->isTightMuon(*primaryVertex));
        eventTuple().other_lepton_muonId.push_back(muonId.GetResultBits());
        if(isMC) {
            const auto match = analysis::gen_truth::LeptonGenMatch(analysis::LorentzVectorM(muon.GetMomentum()), *genParticles);
            eventTuple().other_lepton_gen_match.push_back(static_cast<int>(match.match));
            const auto matched_p4 = match.gen_particle ? match.gen_particle->p4() : analysis::LorentzVectorXYZ();
            eventTuple().other_lepton_gen_p4.push_back(ntuple::LorentzVectorM(matched_p4));
        }
    }
}

void BaseTupleProducer::FillLegGenMatch(const analysis::LorentzVectorXYZ& p4)
{
    using namespace analysis;
    static constexpr int default_int_value = ntuple::DefaultFillValue<Int_t>();

    if(isMC) {
        const auto match = gen_truth::LeptonGenMatch(analysis::LorentzVectorM(p4), *genParticles);
        eventTuple().lep_gen_match.push_back(static_cast<int>(match.match));
        const auto matched_p4 = match.gen_particle ? match.gen_particle->p4() : LorentzVectorXYZ();
        eventTuple().lep_gen_p4.push_back(ntuple::LorentzVectorM(matched_p4));
        const auto matched_visible_p4 = match.visible_daughters_p4;
        eventTuple().lep_gen_visible_p4.push_back(ntuple::LorentzVectorM(matched_visible_p4));
        eventTuple().lep_gen_chargedParticles.push_back(match.n_chargedParticles);
        eventTuple().lep_gen_neutralParticles.push_back(match.n_neutralParticles);
    } else {
        eventTuple().lep_gen_match.push_back(default_int_value);
        eventTuple().lep_gen_p4.push_back(ntuple::LorentzVectorM());
        eventTuple().lep_gen_visible_p4.push_back(ntuple::LorentzVectorM());
        eventTuple().lep_gen_chargedParticles.push_back(default_int_value);
        eventTuple().lep_gen_neutralParticles.push_back(default_int_value);
    }
}

void BaseTupleProducer::FillMetFilters(analysis::Period period)
{
    using MetFilters = ntuple::MetFilters;
    using Filter = MetFilters::Filter;

    MetFilters filters;
    const auto setResult = [&](Filter filter, const std::string& name) {
        bool result;
        auto iter = customMetFilters_token.find(name);
        if(iter != customMetFilters_token.end()) {
            edm::Handle<bool> result_handle;
            edmEvent->getByToken(iter->second, result_handle);
            result = *result_handle;
        }
        else if(!triggerTools.TryGetAnyTriggerResult(name, result))
            throw cms::Exception("TauTriggerSelectionFilter") << "MET filter '" << name << "' not found.";
        filters.SetResult(filter, result);
    };

    setResult(Filter::PrimaryVertex, "Flag_goodVertices");
    setResult(Filter::BeamHalo, "Flag_globalSuperTightHalo2016Filter");
    setResult(Filter::HBHE_noise, "Flag_HBHENoiseFilter");
    setResult(Filter::HBHEiso_noise, "Flag_HBHENoiseIsoFilter");
    setResult(Filter::ECAL_TP, "Flag_EcalDeadCellTriggerPrimitiveFilter");
    setResult(Filter::ee_badSC_noise, "Flag_eeBadScFilter");
    setResult(Filter::badMuon, "Flag_BadPFMuonFilter");

    // if(period == analysis::Period::Run2016){
    //     edm::Handle<bool> badPFMuon;
    //     edmEvent->getByToken(badPFMuonFilter_token, badPFMuon);
    //     filters.SetResult(Filter::badMuon, *badPFMuon);
    //
    //     edm::Handle<bool> badChCandidate;
    //     edmEvent->getByToken(badChCandidateFilter_token, badChCandidate);
    //     filters.SetResult(Filter::badChargedHadron,*badChCandidate);
    // }

    if(period == analysis::Period::Run2017 || period == analysis::Period::Run2018){
        //setResult(Filter::badChargedHadron, "Flag_BadChargedCandidateFilter");
        setResult(Filter::ecalBadCalib, "ecalBadCalibReducedMINIAODFilter");
    }


    eventTuple().metFilters = filters.FilterResults();
}

void BaseTupleProducer::ApplyBaseSelection(analysis::SelectionResultsBase& selection)
{
    using namespace cuts::btag_2016;

    selection.jets = CollectJets();
    if(applyRecoilCorr)
        ApplyRecoilCorrection(selection.jets);

}

std::vector<BaseTupleProducer::ElectronCandidate> BaseTupleProducer::CollectVetoElectrons(
        bool isTightSelection, const std::vector<const ElectronCandidate*>& signalElectrons)
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectVetoElectron, this, _1, _2, signalElectrons, isTightSelection);
    return CollectObjects("vetoElectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::MuonCandidate> BaseTupleProducer::CollectVetoMuons(
        bool isTightSelection, const std::vector<const MuonCandidate*>& signalMuons)
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&BaseTupleProducer::SelectVetoMuon, this, _1, _2, signalMuons, isTightSelection);
    return CollectObjects("vetoMuons", base_selector, muons);
}

std::vector<BaseTupleProducer::JetCandidate> BaseTupleProducer::CollectJets()
{
    using namespace cuts::btag_2016;
    using namespace std::placeholders;
    const auto baseSelector = std::bind(&BaseTupleProducer::SelectJet, this, _1, _2);

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
                                           const std::vector<const ElectronCandidate*>& signalElectrons,
                                           bool isTightSelection) const
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
    const float passID = isTightSelection ? electron->electronID("mvaEleID-Fall17-iso-V2-wp80")
                                          : electron->electronID("mvaEleID-Fall17-iso-V2-wpLoose");
    cut(passID > 0.5, "electronId");
    for(size_t n = 0; n < signalElectrons.size(); ++n) {
        std::ostringstream ss_name;
        ss_name << "isNotSignal_" << n + 1;
        const bool isNotSignal =  &(*electron) != &(*(*signalElectrons.at(n)));
        cut(isNotSignal, ss_name.str());
    }
}

void BaseTupleProducer::SelectVetoMuon(const MuonCandidate& muon, Cutter& cut,
                                       const std::vector<const MuonCandidate*>& signalMuons,
                                       bool isTightSelection) const
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
    double iso_cut = isTightSelection ? tightIso : pfRelIso04;
    cut(muon.GetIsolation() < iso_cut, "iso", muon.GetIsolation());
    bool passMuonId = isTightSelection ? muon->isTightMuon(*primaryVertex) : muon->isLooseMuon();
    cut(passMuonId, "muonID");
    for(size_t n = 0; n < signalMuons.size(); ++n) {
        std::ostringstream ss_name;
        ss_name << "isNotSignal_" << n + 1;
        const bool isNotSignal =  &(*muon) != &(*(*signalMuons.at(n)));
        cut(isNotSignal, ss_name.str());
    }
}

void BaseTupleProducer::SelectJet(const JetCandidate& jet, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::jetID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = jet.GetMomentum();
    cut(p4.Pt() > pt - pt_safety, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < cuts::hh_bbtautau_2017::jetID::eta, "eta", p4.Eta());
    cut(PassPFTightId(*jet,period), "jet_id");
    // for(size_t n = 0; n < signalLeptonMomentums.size(); ++n) {
    //     std::ostringstream cut_name;
    //     cut_name << "deltaR_lep" << n + 1;
    //     const double deltaR = ROOT::Math::VectorUtil::DeltaR(p4, signalLeptonMomentums.at(n));
    //     cut(deltaR > deltaR_signalObjects, cut_name.str(), deltaR);
    // }
}

bool BaseTupleProducer::PassMatchOrIsoSelection(const TauCandidate& tau) const
{
    bool match_condition = false;
    if(isMC){
        const auto match = analysis::gen_truth::LeptonGenMatch(analysis::LorentzVectorM(tau.GetMomentum()), *genParticles);
        match_condition = match.match == analysis::GenLeptonMatch::Tau;
    }
    return (tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ||
        tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5 ||
        tau->tauID("byVVVLooseDeepTau2017v2p1VSjet") > 0.5 ||
        tau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT2016") > 0.5) || match_condition;
}


void BaseTupleProducer::FillElectron(const analysis::SelectionResultsBase& selection)
{
    static const float default_value = ntuple::DefaultFillValue<float>();
    static const float default_value_int = ntuple::DefaultFillValue<int>();
    for(const ElectronCandidate& electron : selection.electrons){
        eventTuple().lep_p4.push_back(ntuple::LorentzVectorM(electron.GetMomentum()));
        eventTuple().lep_q.push_back(electron.GetCharge());
        eventTuple().lep_type.push_back(static_cast<Int_t>(analysis::LegType::e));
        eventTuple().lep_dxy.push_back(electron->gsfTrack()->dxy(primaryVertex->position()));
        eventTuple().lep_dz.push_back(electron->gsfTrack()->dz(primaryVertex->position()));
        eventTuple().lep_iso.push_back(electron.GetIsolation());
        eventTuple().lep_decayMode.push_back(default_value_int);
        eventTuple().lep_oldDecayModeFinding.push_back(false);
        eventTuple().lep_newDecayModeFinding.push_back(false);
        eventTuple().lep_elePassConversionVeto.push_back(electron->passConversionVeto());
        analysis::DiscriminatorIdResults eleId_iso;
        eleId_iso.SetResult(analysis::DiscriminatorWP::Loose,electron->electronID("mvaEleID-Fall17-iso-V2-wpLoose") > 0.5);
        eleId_iso.SetResult(analysis::DiscriminatorWP::Medium,electron->electronID("mvaEleID-Fall17-iso-V2-wp90") > 0.5);
        eleId_iso.SetResult(analysis::DiscriminatorWP::Tight,electron->electronID("mvaEleID-Fall17-iso-V2-wp80") > 0.5);
        eventTuple().lep_eleId_iso.push_back(eleId_iso.GetResultBits());
        analysis::DiscriminatorIdResults eleId_noIso;
        eleId_noIso.SetResult(analysis::DiscriminatorWP::Loose,electron->electronID("mvaEleID-Fall17-noIso-V2-wpLoose") > 0.5);
        eleId_noIso.SetResult(analysis::DiscriminatorWP::Medium,electron->electronID("mvaEleID-Fall17-noIso-V2-wp90") > 0.5);
        eleId_noIso.SetResult(analysis::DiscriminatorWP::Tight,electron->electronID("mvaEleID-Fall17-noIso-V2-wp80") > 0.5);
        eventTuple().lep_eleId_noIso.push_back(eleId_noIso.GetResultBits());
        eventTuple().lep_muonId.push_back(0);
        for(const auto& tau_id_entry : analysis::tau_id::GetTauIdDescriptors()) {
            const auto& desc = tau_id_entry.second;
            desc.FillTuple<ntuple::EventTuple,pat::Tau>(eventTuple, nullptr, default_value);
        }
        FillLegGenMatch(electron->p4());
    }

}

void BaseTupleProducer::FillMuon(const analysis::SelectionResultsBase& selection)
{
    static const float default_value = ntuple::DefaultFillValue<float>();
    static const float default_value_int = ntuple::DefaultFillValue<int>();
    for(const MuonCandidate& muon : selection.muons){
        eventTuple().lep_p4.push_back(ntuple::LorentzVectorM(muon.GetMomentum()));
        eventTuple().lep_q.push_back(muon.GetCharge());
        eventTuple().lep_type.push_back(static_cast<Int_t>(analysis::LegType::mu));
        eventTuple().lep_dxy.push_back(muon->muonBestTrack()->dxy(primaryVertex->position()));
        eventTuple().lep_dz.push_back(muon->muonBestTrack()->dz(primaryVertex->position()));
        eventTuple().lep_iso.push_back(muon.GetIsolation());
        eventTuple().lep_decayMode.push_back(default_value_int);
        eventTuple().lep_oldDecayModeFinding.push_back(false);
        eventTuple().lep_newDecayModeFinding.push_back(false);
        eventTuple().lep_elePassConversionVeto.push_back(false);
        eventTuple().lep_eleId_iso.push_back(0);
        eventTuple().lep_eleId_noIso.push_back(0);
        analysis::DiscriminatorIdResults muonId;
        muonId.SetResult(analysis::DiscriminatorWP::Loose,muon->isLooseMuon());
        muonId.SetResult(analysis::DiscriminatorWP::Medium,muon->isMediumMuon());
        muonId.SetResult(analysis::DiscriminatorWP::Tight,muon->isTightMuon(*primaryVertex));
        eventTuple().lep_muonId.push_back(muonId.GetResultBits());
        for(const auto& tau_id_entry : analysis::tau_id::GetTauIdDescriptors()) {
            const auto& desc = tau_id_entry.second;
            desc.FillTuple<ntuple::EventTuple,pat::Tau>(eventTuple, nullptr, default_value);
        }
        FillLegGenMatch(muon->p4());
    }
}

void BaseTupleProducer::FillTau(const analysis::SelectionResultsBase& selection)
{
    static const float default_value = ntuple::DefaultFillValue<float>();
    for(const TauCandidate& tau : selection.taus) {
        eventTuple().lep_p4.push_back(ntuple::LorentzVectorM(tau.GetMomentum()));
        eventTuple().lep_q.push_back(tau.GetCharge());
        eventTuple().lep_type.push_back(static_cast<Int_t>(analysis::LegType::tau));
        const auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
        eventTuple().lep_dxy.push_back(packedLeadTauCand->dxy());
        eventTuple().lep_dz.push_back(packedLeadTauCand->dz());
        eventTuple().lep_iso.push_back(default_value);
        eventTuple().lep_decayMode.push_back(tau->decayMode());
        bool oldDM = tau->tauID("decayModeFinding") > 0.5f;
        eventTuple().lep_oldDecayModeFinding.push_back(oldDM);
        bool newDM = tau->tauID("decayModeFindingNewDMs") > 0.5f;
        eventTuple().lep_newDecayModeFinding.push_back(newDM);
        eventTuple().lep_elePassConversionVeto.push_back(false);
        eventTuple().lep_eleId_iso.push_back(0);
        eventTuple().lep_eleId_noIso.push_back(0);
        eventTuple().lep_muonId.push_back(0);
        for(const auto& tau_id_entry : analysis::tau_id::GetTauIdDescriptors()) {
            const auto& desc = tau_id_entry.second;
            desc.FillTuple(eventTuple, &*tau, default_value); //or tau->getPtr()
        }
        FillLegGenMatch(tau->p4());
    }
}

void BaseTupleProducer::FillHiggsDaughtersIndexes(const analysis::SelectionResultsBase& selection, size_t shift)
{
    for(unsigned n = 0; n < selection.higgses_pair_indexes.size(); ++n){
        const auto higgs_pair = selection.higgses_pair_indexes.at(n);
        eventTuple().first_daughter_indexes.push_back(higgs_pair.first);
        eventTuple().second_daughter_indexes.push_back(shift + higgs_pair.second);
    }
}


void BaseTupleProducer::FillEventTuple(const analysis::SelectionResultsBase& selection)
{
    using namespace analysis;

    eventTuple().run  = edmEvent->id().run();
    eventTuple().lumi = edmEvent->id().luminosityBlock();
    eventTuple().evt  = edmEvent->id().event();
    eventTuple().genEventType = static_cast<int>(GenEventType::Other);
    eventTuple().genEventWeight = isMC ? genEvt->weight() : 1;

    eventTuple().npv = vertices->size();
    eventTuple().npu = gen_truth::GetNumberOfPileUpInteractions(PUInfo);
    eventTuple().rho = *rho;

    // HTT candidate
    for(const auto& result : selection.svfitResult){
        eventTuple().SVfit_Higgs_index.push_back(result.first);
        eventTuple().SVfit_is_valid.push_back(result.second.has_valid_momentum);
        eventTuple().SVfit_p4.push_back(ntuple::LorentzVectorM(result.second.momentum));
        eventTuple().SVfit_p4_error.push_back(ntuple::LorentzVectorM(result.second.momentum_error));
        eventTuple().SVfit_mt.push_back(result.second.transverseMass);
        eventTuple().SVfit_mt_error.push_back(result.second.transverseMass_error);
        eventTuple().SVfit_unc_source.push_back(static_cast<Int_t>(UncertaintySource::None));
        eventTuple().SVfit_unc_scale.push_back(static_cast<Int_t>(UncertaintyScale::Central));
    }


    // MET
    eventTuple().pfMET_p4 = met->GetMomentum();
    eventTuple().pfMET_cov = met->GetCovMatrix();
    if(isMC & !isEmbedded){
      ntuple::LorentzVectorM genMet_momentum(genMET->at(0).pt(),0,genMET->at(0).eta(),0);
      eventTuple().genMET_p4 = genMet_momentum;
    }

    FillMetFilters(period);

    std::set<const pat::Jet*> selected_jets;

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
        analysis::DiscriminatorIdResults jet_pu_id;
        jet_pu_id.SetResult(analysis::DiscriminatorWP::Loose,jet->userInt("pileupJetId:fullId") & (1 << 2));
        jet_pu_id.SetResult(analysis::DiscriminatorWP::Medium,jet->userInt("pileupJetId:fullId") & (1 << 1));
        jet_pu_id.SetResult(analysis::DiscriminatorWP::Tight,jet->userInt("pileupJetId:fullId") & (1 << 0));
        eventTuple().jets_pu_id.push_back(jet_pu_id.GetResultBits());
        eventTuple().jets_hadronFlavour.push_back(jet->hadronFlavour());
        // Jet resolution
        JME::JetParameters parameters;
        parameters.setJetPt(jet.GetMomentum().pt());
        parameters.setJetEta(jet.GetMomentum().eta());
        parameters.setRho(*rho);
        float jet_resolution = resolution.getResolution(parameters);
        eventTuple().jets_resolution.push_back(jet_resolution); // percentage

        const auto raw_match_bits = triggerTools.GetJetMatchBits(p4, cuts::H_tautau_2016::DeltaR_triggerMatch);
        const auto match_bits = TriggerDescriptorCollection::ConvertToRootRepresentation(raw_match_bits);
        for(size_t n = 0; n < match_bits.size(); ++n) {
            const std::string br_name = "jets_triggerFilterMatch_" + std::to_string(n);
            eventTuple.get<std::vector<ULong64_t>>(br_name).push_back(match_bits.at(n));
        }
    }

    for(const auto jet_cand : jets){
        const auto pat_jet = &(*jet_cand);
        if(selected_jets.count(pat_jet)) continue;
        const LorentzVector& other_p4 = jet_cand.GetMomentum();
        eventTuple().other_jets_p4.push_back(ntuple::LorentzVectorE(other_p4));
    }




    for(const JetCandidate& jet : fatJets) {
        eventTuple().fatJets_p4.push_back(ntuple::LorentzVectorE(jet.GetMomentum()));

        std::string subjets_collection;
        if(period == Period::Run2016) {
            subjets_collection = "SoftDrop";
            eventTuple().fatJets_m_softDrop.push_back(GetUserFloat(jet, "ak8PFJetsCHSSoftDropMass"));
            eventTuple().fatJets_jettiness_tau1.push_back(GetUserFloat(jet, "NjettinessAK8:tau1"));
            eventTuple().fatJets_jettiness_tau2.push_back(GetUserFloat(jet, "NjettinessAK8:tau2"));
            eventTuple().fatJets_jettiness_tau3.push_back(GetUserFloat(jet, "NjettinessAK8:tau3"));
        } else if(period == Period::Run2017 ||period == Period::Run2018 ) {
            subjets_collection = "SoftDropPuppi";
            eventTuple().fatJets_m_softDrop.push_back(GetUserFloat(jet, "ak8PFJetsPuppiSoftDropMass"));
            eventTuple().fatJets_jettiness_tau1.push_back(GetUserFloat(jet, "NjettinessAK8Puppi:tau1"));
            eventTuple().fatJets_jettiness_tau2.push_back(GetUserFloat(jet, "NjettinessAK8Puppi:tau2"));
            eventTuple().fatJets_jettiness_tau3.push_back(GetUserFloat(jet, "NjettinessAK8Puppi:tau3"));
            eventTuple().fatJets_jettiness_tau4.push_back(GetUserFloat(jet, "NjettinessAK8Puppi:tau4"));
        }
        else
            throw analysis::exception("FillEventTuple: fill fatJets for period %1% is not implemented") % period;

        if(!jet->hasSubjets(subjets_collection)) continue;
        const size_t parentIndex = eventTuple().fatJets_p4.size() - 1;
        const auto& sub_jets = jet->subjets(subjets_collection);
        for(const auto& sub_jet : sub_jets) {
            eventTuple().subJets_p4.push_back(ntuple::LorentzVectorE(sub_jet->p4()));
            eventTuple().subJets_parentIndex.push_back(parentIndex);
        }
    }

    FillLheInfo();
    if(isMC) {
        FillGenParticleInfo();
        FillGenJetInfo();
    }

    eventTuple().extraelec_veto = selection.electronVeto;
    eventTuple().extramuon_veto = selection.muonVeto;
    FillOtherLeptons(selection.other_electrons,selection.other_muons);

    // eventTuple().trigger_match = !applyTriggerMatch || selection.triggerResults.AnyAcceptAndMatch();
    eventTuple().trigger_accepts = selection.triggerResults.at(0).GetAcceptBits();
    for(unsigned n = 0; n < selection.triggerResults.size(); ++n){
        eventTuple().trigger_matches.push_back(selection.triggerResults.at(n).GetMatchBits());
    }

}
