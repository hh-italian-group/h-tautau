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
        triggerTools.SetTriggerAcceptBits(selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
    }

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

    if(applyTriggerMatch){
        triggerTools.SetTriggerMatchBits(selection.triggerResults, selected_higgs,
                                         cuts::H_tautau_2016::DeltaR_triggerMatch);
        cut(selection.triggerResults.AnyAcceptAndMatch(), "trigger_match");
    }


    selection.SetHiggsCandidate(selected_higgs);

    //Third-Lepton Veto
    selection.other_electrons = CollectVetoElectrons({ &selection.higgs->GetFirstDaughter() });
    selection.other_muons = CollectVetoMuons();
    selection.electronVeto = selection.other_electrons.size();
    selection.muonVeto = selection.other_muons.size();

    ApplyBaseSelection(selection, selection.higgs->GetDaughterMomentums());
    if(runSVfit)
        selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);
    FillEventTuple(selection);
    if(eventEnergyScale == analysis::EventEnergyScale::Central)
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
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

void TupleProducer_eTau::SelectSignalElectron(const ElectronCandidate& electron, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::electronID;

    cut(true, "gt0_cand");
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
        cut(electron->passConversionVeto(), "conversionVeto");
    } else {
        if(period != analysis::Period::Run2017)
            cut(electron.GetIsolation() < pfRelIso04, "iso", electron.GetIsolation());
    }
}

void TupleProducer_eTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::tauID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    double pt_cut = pt;
    if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::ETau::tauID::pt;
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
