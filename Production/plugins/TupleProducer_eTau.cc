/*! Implementation of an event tuple producer for the e-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_eTau.h"
#include "../interface/GenTruthTools.h"
#include "h-tautau/Cuts/include/hh_bbtautau_Run2.h"

void TupleProducer_eTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::hh_bbtautau_Run2::ETau;
    using HiggsCandidate = analysis::CompositeCandidate<ElectronCandidate,TauCandidate>;

    SelectionResultsBase selection(eventId);
    cut(primaryVertex.isNonnull(), "vertex");

    analysis::TriggerResults refTriggerResults;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        if(applyTriggerMatchCut) cut(refTriggerResults.AnyAccept(), "trigger");
    }

    // Signal-like leptons selection
    auto electrons = CollectSignalElectrons();
    cut(electrons.size(), "electrons");

    std::sort(electrons.begin(),electrons.end(),&analysis::LeptonComparitor<ElectronCandidate>);
    selection.electrons.push_back(electrons.at(0));
    selection.other_electrons = CollectVetoElectrons(false,{&electrons.at(0)});
    selection.electronVeto = selection.other_electrons.size();

    selection.other_muons = CollectVetoMuons(false,{});
    selection.muonVeto = selection.other_muons.size();

    auto other_tight_electrons = CollectVetoElectrons(true,{&electrons.at(0)});
    cut(other_tight_electrons.empty(), "tightElectronVeto");

    auto other_tight_muons = CollectVetoMuons(true,{});
    cut( other_tight_muons.empty(), "tightElectronVeto");

    selection.taus = CollectSignalTaus();
    cut(selection.taus.size(), "taus");

    static constexpr double DeltaR_betweenSignalObjects = cuts::hh_bbtautau_Run2::DeltaR_Lep_Lep;
    auto higgses_indexes = FindCompatibleObjects(selection.electrons, selection.taus, DeltaR_betweenSignalObjects,
                                                 "H_e_tau");
    cut(higgses_indexes.size(), "ele_tau_pair");

    for(size_t n = 0; n < higgses_indexes.size(); ++n){
        auto daughter_index = higgses_indexes.at(n);
        const HiggsCandidate selected_higgs(selection.electrons.at(daughter_index.first),
                                            selection.taus.at(daughter_index.second));

        if(applyTriggerMatch){
            analysis::TriggerResults triggerResults(refTriggerResults);
            triggerTools.SetTriggerMatchBits(triggerResults, selected_higgs,
                                             cuts::hh_bbtautau_Run2::DeltaR_triggerMatch);
            selection.triggerResults.push_back(triggerResults);
        }

        selection.higgses_pair_indexes.push_back(daughter_index);
    }

    ApplyBaseSelection(selection);
    FillEventTuple(selection);
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
    using namespace cuts::hh_bbtautau_Run2::ETau::electronID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = electron.GetMomentum();
    cut(p4.pt() > pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double electron_xy = std::abs(electron->gsfTrack()->dxy(primaryVertex->position()));
    cut(electron_xy < dxy, "dxy", electron_xy);
    const double electron_dz = std::abs(electron->gsfTrack()->dz(primaryVertex->position()));
    cut(electron_dz < dz, "dz", electron_dz);
    const bool passID = electron->electronID(cuts::electronID_Run2::mvaEleID_iso_Tight) > 0.5f ||
                        electron->electronID(cuts::electronID_Run2::mvaEleID_noIso_Tight) > 0.5f;
    cut(passID, "electronId");
}

void TupleProducer_eTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    cut(true, "gt0_cand");
    cut(PassMatchSelection(tau) || PassIsoSelection(tau), "iso");
}

void TupleProducer_eTau::FillEventTuple(const SelectionResultsBase& selection)
{
    using Channel = analysis::Channel;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;

    Lock lock(eventTuple.GetMutex());

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelId = static_cast<int>(Channel::ETau);

    BaseTupleProducer::FillElectron(selection);
    BaseTupleProducer::FillTau(selection);
    BaseTupleProducer::FillHiggsDaughtersIndexes(selection,selection.electrons.size());

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_eTau);
