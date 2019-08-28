/*! Implementation of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_muTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_muTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::MuTau;

    SelectionResultsBase selection(eventId);
    cut(primaryVertex.isNonnull(), "vertex");

    analysis::TriggerResults refTriggerResults;
    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        if(applyTriggerMatchCut) cut(refTriggerResults.AnyAccept(), "trigger");
    }

    // Signal-like leptons selection
    auto muons = CollectSignalMuons();
    cut(muons.size(), "muons");

    std::sort(muons.begin(),muons.end(),&LeptonComparitor<MuonCandidate>);
    selection.muons.push_back(muons.at(0));
    selection.other_muons = CollectVetoMuons({&muons.at(0)});
    selection.muonVeto = selection.other_muons.size();

    selection.other_electrons = CollectVetoElectrons();
    selection.electronVeto = selection.other_electrons.size();

    auto other_tight_electrons = CollectVetoElectrons(true);
    cut(other_tight_electrons.empty(), "tightElectronVeto");

    auto other_tight_muons = CollectVetoMuons(true,{&muons.at(0)});
    cut(other_tight_muons.empty(), "tightElectronVeto");

    selection.taus = CollectSignalTaus();
    cut(selection.taus.size(), "taus");

    static constexpr double DeltaR_betweenSignalObjects = cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects;

    auto higgses_indexes = FindCompatibleObjects(selection.muons, selection.taus, DeltaR_betweenSignalObjects, "H_mu_tau");
    cut(higgses_indexes.size(), "mu_tau_pair");

    for(size_t n = 0; n < higgses_indexes.size(); ++n){
        auto daughter_index = higgses_indexes.at(n);
        analysis::CompositeCandidate<MuonCandidate,TauCandidate> selected_higgs =
            analysis::CompositeCandidate<MuonCandidate,TauCandidate>(selection.muons.at(daughter_index.first), selection.taus.at(daughter_index.second));

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

}

std::vector<BaseTupleProducer::MuonCandidate> TupleProducer_muTau::CollectSignalMuons()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muTau::SelectSignalMuon, this, _1, _2);
    return CollectObjects("SignalMuons", base_selector, muons);
}

std::vector<BaseTupleProducer::TauCandidate> TupleProducer_muTau::CollectSignalTaus()
{
    using namespace std::placeholders;
    const auto base_selector = std::bind(&TupleProducer_muTau::SelectSignalTau, this, _1, _2);
    return CollectObjects("SignalTaus", base_selector, taus);
}

void TupleProducer_muTau::SelectSignalMuon(const MuonCandidate& muon, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::muonID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = muon.GetMomentum();
    static constexpr double pt_cut = cuts::hh_bbtautau_2017::MuTau::muonID::pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    cut(std::abs(p4.eta()) < eta, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    cut(muon->isMediumMuon(), "muonID");
}

void TupleProducer_muTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::tauID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    cut(p4.Pt() > pt - BaseTupleProducer::pt_shift , "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    bool match_condition = false;
    if(isMC){
        const auto match = analysis::gen_truth::LeptonGenMatch(analysis::LorentzVectorM(p4), *genParticles);
        match_condition = match.match == analysis::GenLeptonMatch::Tau;
    }
    bool iso_condition = (tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 ||
        tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5 ||
        tau->tauID("byVVVLooseDeepTau2017v2VSjet") > 0.5) || match_condition;
    cut(iso_condition, "iso");
}


void TupleProducer_muTau::FillEventTuple(const SelectionResultsBase& selection)
{
    using Channel = analysis::Channel;
    using Mutex = std::recursive_mutex;
    using Lock = std::lock_guard<Mutex>;

    Lock lock(eventTuple.GetMutex());

    BaseTupleProducer::FillEventTuple(selection);
    eventTuple().channelId = static_cast<int>(Channel::MuTau);

    BaseTupleProducer::FillMuon(selection);
    BaseTupleProducer::FillTau(selection);
    BaseTupleProducer::FillHiggsDaughtersIndexes(selection,selection.muons.size());

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muTau);
