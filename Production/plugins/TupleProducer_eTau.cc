/*! Implementation of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_eTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_eTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::ETau;

    SelectionResultsBase selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    analysis::TriggerResults refTriggerResults;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        if(applyTriggerMatchCut) cut(refTriggerResults.AnyAccept(), "trigger");
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

    const double DeltaR_betweenSignalObjects = cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects;
    auto higgses_indexes = FindCompatibleObjects(selection.electrons, selection.taus, DeltaR_betweenSignalObjects, "H_e_tau");
    cut(higgses_indexes.size(), "ele_tau_pair");

    for(size_t n = 0; n < higgses_indexes.size(); ++n){
        auto daughter_index = higgses_indexes.at(n);
        analysis::CompositeCandidate<ElectronCandidate,TauCandidate> selected_higgs =
            analysis::CompositeCandidate<ElectronCandidate,TauCandidate>(selection.electrons.at(daughter_index.first), selection.taus.at(daughter_index.second));

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
        previous_selection = SelectionResultsBasePtr(new SelectionResultsBase(selection));
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
    double pt_cut = cuts::hh_bbtautau_2017::ETau::electronID::pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electron_xy = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electron_xy < dxy, "dxy", electron_xy);
    const double electron_dz = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electron_dz < dz, "dz", electron_dz);
    float passID = electron->electronID("mvaEleID-Fall17-noIso-V1-wp90");
    cut(passID == 1, "electronId");
}

void TupleProducer_eTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::ETau::tauID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    double pt_cut = pt;
    cut(p4.Pt() > pt_cut - BaseTupleProducer::pt_shift, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
}

void TupleProducer_eTau::FillEventTuple(const SelectionResultsBase& selection)
{
    using Channel = analysis::Channel;
  
    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::ETau);

    BaseTupleProducer::FillElectron(selection);
    BaseTupleProducer::FillTau(selection);
    BaseTupleProducer::FillHiggsDaughtersIndexes(selection,selection.electrons.size());

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_eTau);
