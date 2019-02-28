/*! Implementation of an event tuple producer for the mu-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_muTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_muTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::MuTau;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) {
        TriggerResults refTriggerResults;
        triggerTools.SetTriggerAcceptBits(refTriggerResults);
        cut(refTriggerResults.AnyAccpet(), "trigger");
    }

    // Signal-like leptons selection
    auto muons = CollectSignalMuons();
    cut(muons.size(), "muons");

    std::sort(muons.begin(),muons.end(),&LeptonComparitor<MuonCandidate>);
    selection.muons.push_back(muons.at(0));
    selection.other_muons = CollectVetoMuons({&muons.at(0)});
    selection.muonVeto = selection.other_muons.size();
    cut(!selection.muonVeto, "no_extra_muon");

    selection.other_electrons = CollectVetoElectrons();
    selection.electronVeto = selection.other_electrons.size();
    cut(!selection.electronVeto, "no_extra_ele");

    selection.taus = CollectSignalTaus();
    cut(selection.taus.size() > 0, "taus");



    const double DeltaR_betweenSignalObjects = (productionMode == ProductionMode::hh ||
            productionMode == ProductionMode::tau_pog)
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selection.muons, selection.taus, DeltaR_betweenSignalObjects, "H_mu_tau");
    cut(higgses.size(), "mu_tau_pair");

    for(size_t n = 0; n < higgses.size(); ++n){
        HiggsCandidate selected_higgs = higgses.at(n);
        if (higgses.at(n).GetFirstDaughter().GetMomentum().Pt() < higgses.at(n).GetSecondDaughter().GetMomentum().Pt())
            selected_higgs = HiggsCandidate(higgses.at(n).GetSecondDaughter(), higgses.at(n).GetFirstDaughter());

        if(applyTriggerMatch){
            TriggerResults triggerResults(refTriggerResults);
            triggerTools.SetTriggerMatchBits(triggerResults, selected_higgs,
                                          cuts::H_tautau_2016::DeltaR_triggerMatch);
            selection.triggerResults.push_back(triggerResults);
        }
        selection.SetHiggsCandidate(selected_higgs);
        selection.higgses.push_back(selected_higgs);

        if(runSVfit)
            selection.svfitResult = svfitProducer->Fit(*selection.higgs, *met);

    }

    ApplyBaseSelection(selection);

    FillEventTuple(selection);

    if(eventEnergyScale == analysis::EventEnergyScale::Central)
        previous_selection = SelectionResultsPtr(new SelectionResults(selection));
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
    double pt_cut = pt;
    if(productionMode == ProductionMode::hh || productionMode == ProductionMode::tau_pog) {
        pt_cut = period == analysis::Period::Run2017 ? cuts::hh_bbtautau_2017::MuTau::muonID::pt : cuts::hh_bbtautau_2016::MuTau::muonID::pt;
    }
    else if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::MuTau::muonID::pt;
    else if (productionMode == ProductionMode::h_tt_sm) pt_cut = cuts::H_tautau_2016_sm::MuTau::muonID::pt;
    cut(p4.pt() > pt_cut, "pt", p4.pt());
    const double eta_cut  = productionMode == ProductionMode::h_tt_sm ? cuts::H_tautau_2016_sm::MuTau::muonID::eta : eta;
    cut(std::abs(p4.eta()) < eta_cut, "eta", p4.eta());
    const double muon_dxy = std::abs(muon->muonBestTrack()->dxy(primaryVertex->position()));
    cut(muon_dxy < dxy, "dxy", muon_dxy);
    const double muon_dz = std::abs(muon->muonBestTrack()->dz(primaryVertex->position()));
    cut(muon_dz < dz, "dz", muon_dz);
    if(productionMode == ProductionMode::hh){
        cut(muon->isTightMuon(*primaryVertex), "muonID");
        cut(muon.GetIsolation() < pfRelIso04, "iso", muon.GetIsolation());
    }
    else cut(muon->isMediumMuon(), "muonID");

}

void TupleProducer_muTau::SelectSignalTau(const TauCandidate& tau, Cutter& cut) const
{
    using namespace cuts::H_tautau_2016::MuTau::tauID;

    cut(true, "gt0_cand");
    const LorentzVector& p4 = tau.GetMomentum();
    double pt_cut = pt;
    if (productionMode == ProductionMode::h_tt_mssm) pt_cut = cuts::H_tautau_2016_mssm::MuTau::tauID::pt;
    cut(p4.Pt() > pt_cut, "pt", p4.Pt());
    cut(std::abs(p4.Eta()) < eta, "eta", p4.Eta());
    if(productionMode == ProductionMode::hh){
        const auto dmFinding = tau->tauID("decayModeFinding");
        cut(dmFinding > decayModeFinding, "decayMode", dmFinding);
    }
    auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronVLooseMVA6") > againstElectronVLooseMVA6, "againstElectron");
        cut(tau->tauID("againstMuonTight3") > againstMuonTight3, "againstMuon");
        if(period == analysis::Period::Run2017) {
            cut(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5, "VVLooseIso");
        }
    }
}

void TupleProducer_muTau::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::MuTau);

    ntuple::StorageMode storageMode(eventTuple().storageMode);
    eventTuple().storageMode = storageMode.Mode();

    BaseTupleProducer::FillMuon(selection);
    BaseTupleProducer::FillTau(selection);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_muTau);
