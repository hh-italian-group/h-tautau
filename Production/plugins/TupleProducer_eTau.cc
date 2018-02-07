/*! Implementation of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_eTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_eTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::ETau;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(triggerDescriptors, selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
    }

    //Di-Lepton Veto
    const auto z_electrons = CollectZelectrons();
    const auto z_electrons_candidates = FindCompatibleObjects(z_electrons, z_electrons, ZeeVeto::deltaR, "Z_e_e",0);
    selection.Zveto = z_electrons_candidates.size();

    // Signal-like leptons selection
    const auto selectedElectrons = CollectSignalElectrons();
    cut(selectedElectrons.size(), "electrons");
    const auto selectedTaus = CollectSignalTaus();
    cut(selectedTaus.size(), "taus");

    const double DeltaR_betweenSignalObjects = productionMode == ProductionMode::hh
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selectedElectrons, selectedTaus, DeltaR_betweenSignalObjects, "H_e_tau");
    cut(higgses.size(), "ele_tau_pair");

    std::sort(higgses.begin(), higgses.end(), &HiggsComparitor<HiggsCandidate>);
    auto selected_higgs = higgses.front();

    if(applyTriggerMatch)
        triggerTools.SetTriggerMatchBits(triggerDescriptors, selection.triggerResults, selected_higgs,
                                         cuts::H_tautau_2016::DeltaR_triggerMatch, false);

    selection.SetHiggsCandidate(selected_higgs);

    //Third-Lepton Veto
    const auto electronVetoCollection = CollectVetoElectrons({ &selection.higgs->GetFirstDaughter() });
    const auto muonVetoCollection = CollectVetoMuons();
    selection.electronVeto = electronVetoCollection.size();
    selection.muonVeto = muonVetoCollection.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);
    FillEventTuple(selection);
    if(eventEnergyScale == analysis::EventEnergyScale::Central)
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
}

bool TupleProducer_eTau::SelectSpring15VetoElectron(const pat::Electron& electron)
{
    double full5x5_sigmaIetaIeta = electron.full5x5_sigmaIetaIeta();
    double dEtaIn = std::abs(electron.deltaEtaSuperClusterTrackAtVtx());
    float dPhiIn = electron.deltaPhiSuperClusterTrackAtVtx();
    float hOverE = electron.hadronicOverEm();
    const float ecal_energy_inverse = 1.0/electron.ecalEnergy();
    const float eSCoverP = electron.eSuperClusterOverP();
    float ooEmooP =  std::abs(1.0 - eSCoverP)*ecal_energy_inverse;
    constexpr reco::HitPattern::HitCategory missingHitType =
    reco::HitPattern::MISSING_INNER_HITS;
    const unsigned mHits = 
    electron.gsfTrack()->hitPattern().numberOfAllHits(missingHitType);
    bool result = false;
    if(fabs(electron.superCluster()->position().eta()) <= 1.479){
	result=full5x5_sigmaIetaIeta < 0.0114 &&
               fabs(dEtaIn) < 0.0152  &&
	       fabs(dPhiIn) < 0.216   &&
	       hOverE < 0.181 &&
               ooEmooP < 0.207 &&
               mHits <= 2 &&
               electron.passConversionVeto();
    }
    else if(fabs(electron.superCluster()->position().eta()) > 1.479 && fabs(electron.superCluster()->position().eta()) < 2.5){
	result=full5x5_sigmaIetaIeta < 0.0352 &&
               fabs(dEtaIn) < 0.0113  &&
	       fabs(dPhiIn) < 0.237   &&
	       hOverE < 0.116 &&
               ooEmooP < 0.174 &&
               mHits <= 3 &&
               electron.passConversionVeto();       
    } 
    return result;
}

std::vector<BaseTupleProducer::ElectronCandidate> TupleProducer_eTau::CollectZelectrons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_eTau::SelectZElectron, this, _1, _2);
    return CollectObjects("Zelectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::ElectronCandidate> TupleProducer_eTau::CollectSignalElectrons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_eTau::SelectSignalElectron, this, _1, _2);
    return CollectObjects("SignalElectrons", base_selector, electrons);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_eTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_eTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_eTau::SelectZElectron(const ElectronCandidate& electron, Cutter& cut) const
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
    const bool veto  = (*loose_id_veto)[electron.getPtr()];
//    const bool veto  = SelectSpring15VetoElectron(*electron);
    cut(veto, "cut_based_veto");
    cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
}


void TupleProducer_eTau::SelectSignalElectron(const ElectronCandidate& electron, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::electronID;

    cut(true, "gt0_ele_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    double pt_cut = pt;
    if( productionMode == ProductionMode::hh) {
        pt_cut = period == analysis::Period::Run2017 ? cuts::hh_bbtautau_2017::ETau::electronID::pt : cuts::hh_bbtautau_2016::ETau::electronID::pt;
    }
    else if(productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::ETau::electronID::pt;
    else if(productionMode == ProductionMode::h_tt_sm) pt_cut = cuts::H_tautau_2016_sm::ETau::electronID::pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electron_xy = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electron_xy < dxy, "dxy", electron_xy);
    const double electron_dz = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electron_dz < dz, "dz", electron_dz);
    const bool isTight = (*tight_id_decisions)[electron.getPtr()];
    cut(isTight, "electronMVATightID");
    if(productionMode != ProductionMode::hh) {
        const auto eleMissingHits =
                electron->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
        cut(eleMissingHits <= missingHits, "missingHits", eleMissingHits);
        cut(electron->passConversionVeto(), "conversionVeto");
    } else {
        cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
    }
}

void TupleProducer_eTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::tauID;

    cut(true, "gt0_tau_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    double pt_cut = pt;
    if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::ETau::tauID::pt;
    if (period == analysis::Period::Run2017) pt_cut = cuts::hh_bbtautau_2017::ETau::tauID::pt;
    cut(p4.Pt() > pt_cut, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    const auto dmFinding = tau->tauID("decayModeFinding");
    cut(dmFinding > decayModeFinding, "decayMode", dmFinding);
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronTightMVA6") > againstElectronTightMVA6, "againstElectron");
        cut(tau->tauID("againstMuonLoose3") > againstMuonLoose3, "againstMuon");
    }
}


void TupleProducer_eTau::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::ETau);

    ntuple::StorageMode storageMode(eventTuple().storageMode);
    const bool store_tauIds = !previous_selection || !selection.HaveSameSecondLegOrigin(*previous_selection);
    storageMode.SetPresence(EventPart::SecondTauIds, store_tauIds);
    eventTuple().storageMode = storageMode.Mode();

    FillElectronLeg(1, selection.higgs->GetFirstDaughter());
    FillTauLeg(2, selection.higgs->GetSecondDaughter(), store_tauIds);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_eTau);
