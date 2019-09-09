/*! Implementation of an event tuple producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_tauTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_tauTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::TauTau;

    SelectionResultsBase selection(eventId);
    cut(primaryVertex.isNonnull(), "vertex");

    analysis::TriggerResults refTriggerResults;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        if(applyTriggerMatchCut) cut(refTriggerResults.AnyAccept(), "trigger");
    }

    //Third-Lepton Veto
    selection.other_electrons = CollectVetoElectrons(false);
    selection.other_muons = CollectVetoMuons(false);
    selection.electronVeto = selection.other_electrons.size();
    selection.muonVeto = selection.other_muons.size();

    auto other_tight_electrons = CollectVetoElectrons(true);
    cut(other_tight_electrons.empty(), "tightElectronVeto");

    auto other_tight_muons = CollectVetoMuons(true);
    cut(other_tight_muons.empty(), "tightElectronVeto");

    selection.taus = CollectSignalTaus();
    cut(selection.taus.size() > 1, "taus");

    static constexpr double DeltaR_betweenSignalObjects = cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects;
    auto higgses_indexes = FindCompatibleObjects(selection.taus, selection.taus, DeltaR_betweenSignalObjects, "H_tau_tau");
    cut(higgses_indexes.size(), "tau_tau_pair");

    for(size_t n = 0; n < higgses_indexes.size(); ++n){
        auto daughter_index = higgses_indexes.at(n);
        analysis::CompositeCandidate<TauCandidate,TauCandidate> selected_higgs = analysis::CompositeCandidate<TauCandidate,TauCandidate>(selection.taus.at(daughter_index.first), selection.taus.at(daughter_index.second));

        if(applyTriggerMatch){
            analysis::TriggerResults triggerResults(refTriggerResults);
            triggerTools.SetTriggerMatchBits(triggerResults, selected_higgs,
                                          cuts::H_tautau_2016::DeltaR_triggerMatch);
            selection.triggerResults.push_back(triggerResults);
        }

        selection.higgses_pair_indexes.push_back(daughter_index);

        if(runSVfit)
            selection.svfitResult[n] = svfitProducer->Fit(selected_higgs, *met);

    }

    ApplyBaseSelection(selection);

    FillEventTuple(selection);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_tauTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_tauTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_tauTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::TauTau::tauID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    static constexpr double pt_cut = cuts::hh_bbtautau_2017::TauTau::tauID::pt;
    cut(p4.Pt() > pt_cut - BaseTupleProducer::pt_shift, "pt", p4.Pt());
    static constexpr double eta_cut = cuts::hh_bbtautau_2017::TauTau::tauID::eta;
    cut(std::abs(p4.Eta()) < eta_cut, "eta", p4.Eta());
    const auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    bool iso_condition = PassMatchOrIsoSelection(tau);
    cut(iso_condition, "iso");
}

void TupleProducer_tauTau::FillEventTuple(const SelectionResultsBase& selection)
{
    using Channel = analysis::Channel;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;

    Lock lock(eventTuple.GetMutex());

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelId = static_cast<int>(Channel::TauTau);

    BaseTupleProducer::FillTau(selection);
    BaseTupleProducer::FillHiggsDaughtersIndexes(selection,0);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_tauTau);
