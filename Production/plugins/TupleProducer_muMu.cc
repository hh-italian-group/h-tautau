/*! Implementation of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_muMu.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_muMu::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::MuMu;

    SelectionResultsBase selection(eventId);
    cut(primaryVertex.isNonnull(), "vertex");

    analysis::TriggerResults refTriggerResults;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        if(applyTriggerMatchCut) cut(refTriggerResults.AnyAccept(), "trigger");
    }

    // Signal-like leptons selection
    auto muons = CollectSignalMuons();
    cut(muons.size() > 1, "muons");

    selection.other_electrons = CollectVetoElectrons();
    selection.electronVeto = selection.other_electrons.size();

    auto other_tight_electrons = CollectVetoElectrons(true);
    cut(other_tight_electrons.empty(), "tightElectronVeto");

    static constexpr double DeltaR_betweenSignalObjects = cuts::hh_bbtautau_2016::MuMu::DeltaR_betweenSignalObjects;
    auto higgses_indexes = FindCompatibleObjects(muons, muons, DeltaR_betweenSignalObjects, "H_mu_mu");
    cut(higgses_indexes.size(), "mu_mu_pair");

    const auto Comparitor = [&](const std::pair<size_t,size_t>& h1, const std::pair<size_t,size_t>& h2) -> auto
    {
        const auto& h1_leg1 = muons.at(h1.first);
        const auto& h2_leg1 = muons.at(h2.first);
        if(h1_leg1 != h2_leg1) {
            if(h1_leg1.GetIsolation() != h2_leg1.GetIsolation()) return h1_leg1.IsMoreIsolated(h2_leg1);
            if(h1_leg1.GetMomentum().pt() != h2_leg1.GetMomentum().pt())
                return h1_leg1.GetMomentum().pt() > h2_leg1.GetMomentum().pt();
        }

        const auto& h1_leg2 = muons.at(h1.second);
        const auto& h2_leg2 = muons.at(h2.second);
        if(h1_leg2 != h2_leg2) {
            if(h1_leg2.GetIsolation() != h2_leg2.GetIsolation()) return h1_leg2.IsMoreIsolated(h2_leg2);
            if(h1_leg2.GetMomentum().pt() != h2_leg2.GetMomentum().pt())
                return h1_leg2.GetMomentum().pt() > h2_leg2.GetMomentum().pt();
        }

        if(h1_leg1.GetMomentum().energy() != h2_leg1.GetMomentum().energy())
            return h1_leg1.GetMomentum().energy() > h2_leg1.GetMomentum().energy();
        if(h1_leg2.GetMomentum().energy() != h2_leg2.GetMomentum().energy())
            return h1_leg2.GetMomentum().energy() > h2_leg2.GetMomentum().energy();

        if(h1_leg1 == h2_leg1 && h1_leg2 == h2_leg2) return false;
        throw analysis::exception("not found a good criteria for best tau pair");
    };

    std::sort(higgses_indexes.begin(), higgses_indexes.end(), Comparitor);
    auto selected_higgs_index = higgses_indexes.front();
    analysis::CompositeCandidate<MuonCandidate,MuonCandidate> selected_higgs(muons.at(selected_higgs_index.first), muons.at(selected_higgs_index.second));

    cut(selected_higgs.GetFirstDaughter().GetIsolation() < muonID::pfRelIso04, "iso_of_1st_daughter");

    selection.muons.push_back(selected_higgs.GetFirstDaughter());
    selection.muons.push_back(selected_higgs.GetSecondDaughter());

    if(applyTriggerMatch){
        analysis::TriggerResults triggerResults(refTriggerResults);
        triggerTools.SetTriggerMatchBits(triggerResults, selected_higgs,
                                      cuts::H_tautau_2016::DeltaR_triggerMatch);
        selection.triggerResults.push_back(triggerResults);
    }

    selection.higgses_pair_indexes.push_back(selected_higgs_index);

    selection.other_muons = CollectVetoMuons(false,{ &selected_higgs.GetFirstDaughter(),
        &selected_higgs.GetSecondDaughter() });
    selection.muonVeto = selection.other_muons.size();

    auto other_tight_muons = CollectVetoMuons(true,{ &selected_higgs.GetFirstDaughter(),
        &selected_higgs.GetSecondDaughter() });
    cut(other_tight_muons.empty(), "no_extra_muon");

    if(runSVfit)
        selection.svfitResult[0] = svfitProducer->Fit(selected_higgs, *met);

    ApplyBaseSelection(selection);

    FillEventTuple(selection);
}

std::vector<BaseTupleProducer::MuonCandidate> TupleProducer_muMu::CollectSignalMuons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muMu::SelectSignalMuon, this, _1, _2);
    return CollectObjects("SignalMuons", base_selector, muons);
}

void TupleProducer_muMu::SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuMu::muonID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    cut(p4.pt() > cuts::hh_bbtautau_2016::MuMu::muonID::pt, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    bool passMuonId = muon->isMediumMuon();
    cut(passMuonId, "muonID");
}

void TupleProducer_muMu::FillEventTuple(const SelectionResultsBase& selection)
{
    using Channel = analysis::Channel;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;

    Lock lock(eventTuple.GetMutex());

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelId = static_cast<int>(Channel::MuMu);

    BaseTupleProducer::FillMuon(selection);
    BaseTupleProducer::FillHiggsDaughtersIndexes(selection,0);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muMu);
