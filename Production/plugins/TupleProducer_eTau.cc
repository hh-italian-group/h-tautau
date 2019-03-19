/*! Implementation of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_eTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_eTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::ETau;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    analysis::TriggerResults refTriggerResults;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        cut(refTriggerResults.AnyAccpet(), "trigger");
    }

    // Signal-like leptons selection
    auto electrons = CollectSignalElectrons();
    cut(electrons.size(), "electrons");

    std::sort(electrons.begin(),electrons.end(),&LeptonComparitor<ElectronCandidate>);
    selection.electrons.push_back(electrons.at(0));
    selection.other_electrons = CollectVetoElectrons({&electrons.at(0)});
    selection.electronVeto = selection.other_electrons.size();
    cut(!selection.electronVeto, "no_extra_electron");

    selection.other_muons = CollectVetoMuons();
    selection.muonVeto = selection.other_muons.size();
    cut(!selection.muonVeto, "no_extra_muon");

    selection.taus = CollectSignalTaus();
    cut(selection.taus.size(), "taus");

    const double DeltaR_betweenSignalObjects = (productionMode == ProductionMode::hh ||
                productionMode == ProductionMode::tau_pog)
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses_indexes = FindCompatibleObjects(selection.electrons, selection.taus, DeltaR_betweenSignalObjects, "H_e_tau");
    cut(higgses_indexes.size(), "ele_tau_pair");

    for(size_t n = 0; n < higgses_indexes.size(); ++n){
        auto daughter_index = higgses_indexes.at(n);
        HiggsCandidate selected_higgs = HiggsCandidate(selection.electrons.at(daughter_index.first), selection.taus.at(daughter_index.second));

        if(applyTriggerMatch){
            analysis::TriggerResults triggerResults(refTriggerResults);
            triggerTools.SetTriggerMatchBits(triggerResults, selected_higgs,
                                          cuts::H_tautau_2016::DeltaR_triggerMatch);
            selection.triggerResults.push_back(triggerResults);
        }

        selection.higgses_pair_indexes.push_back(daughter_index);

        if(runSVfit)
            selection.svfitResult.push_back(svfitProducer->Fit(selected_higgs, *met));

    }

    ApplyBaseSelection(selection);
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
    if( productionMode == ProductionMode::hh || productionMode == ProductionMode::tau_pog) {
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
    if(productionMode == ProductionMode::hh) {
        const auto dmFinding = tau->tauID("decayModeFinding");
        cut(dmFinding > decayModeFinding, "decayMode", dmFinding);
    }
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronTightMVA6") > againstElectronTightMVA6, "againstElectron");
        cut(tau->tauID("againstMuonLoose3") > againstMuonLoose3, "againstMuon");
        if(period == analysis::Period::Run2017) {
            cut(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5, "VVLooseIso");
        }
    }
}

void TupleProducer_eTau::FillHiggsDaughtersIndexes(const SelectionResults& selection)
{
    for(unsigned n = 0; n < selection.higgses_pair_indexes.size(); ++n){
        const auto higgs_pair = selection.higgses_pair_indexes.at(n);
        eventTuple().first_daughter_indexes.push_back(higgs_pair.first);
        eventTuple().second_daughter_indexes.push_back(selection.electrons.size() + higgs_pair.second);
    }
}


void TupleProducer_eTau::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::ETau);

    ntuple::StorageMode storageMode(eventTuple().storageMode);
    eventTuple().storageMode = storageMode.Mode();

    BaseTupleProducer::FillElectron(selection);
    BaseTupleProducer::FillTau(selection);
    FillHiggsDaughtersIndexes(selection);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_eTau);
