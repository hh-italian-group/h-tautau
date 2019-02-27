/*! Implementation of an event tuple producer for the tau-tau channel.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "../interface/TupleProducer_tauTau.h"
#include "../interface/GenTruthTools.h"

void TupleProducer_tauTau::ProcessEvent(Cutter& cut)
{
    using namespace cuts::H_tautau_2016::TauTau;

    SelectionResults selection(eventId, eventEnergyScale);
    cut(primaryVertex.isNonnull(), "vertex");

    if(applyTriggerMatch) {
        triggerTools.SetTriggerAcceptBits(selection.triggerResults);
        cut(selection.triggerResults.AnyAccpet(), "trigger");
    }

    //Third-Lepton Veto
    selection.other_electrons = CollectVetoElectrons();
    selection.other_muons = CollectVetoMuons();
    selection.electronVeto = selection.other_electrons.size();
    selection.muonVeto = selection.other_muons.size();

    cut(!selection.electronVeto, "no_extra_ele");
    cut(!selection.muonVeto, "no_extra_muon");

    selection.taus = CollectSignalTaus();
    cut(selection.taus.size() > 1, "taus");

    const double DeltaR_betweenSignalObjects = (productionMode == ProductionMode::hh ||
        productionMode == ProductionMode::tau_pog)
            ? cuts::hh_bbtautau_2016::DeltaR_betweenSignalObjects
            : cuts::H_tautau_2016::DeltaR_betweenSignalObjects;
    auto higgses = FindCompatibleObjects(selection.taus, selection.taus, DeltaR_betweenSignalObjects, "H_tau_tau");
    cut(higgses.size(), "tau_tau_pair");

    for(size_t n = 0; n < higgses.size(); ++n){
        HiggsCandidate selected_higgs = higgses.at(n);
        if (higgses.at(n).GetFirstDaughter().GetMomentum().Pt() < higgses.at(n).GetSecondDaughter().GetMomentum().Pt())
            selected_higgs = HiggsCandidate(higgses.at(n).GetSecondDaughter(), higgses.at(n).GetFirstDaughter());

        if(applyTriggerMatch){
            triggerTools.SetTriggerMatchBits(selection.triggerResults, selected_higgs,
                                          cuts::H_tautau_2016::DeltaR_triggerMatch);
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
    double pt_cut = period == analysis::Period::Run2017 ? cuts::hh_bbtautau_2017::TauTau::tauID::pt : cuts::H_tautau_2016::TauTau::tauID::pt;
    cut(p4.Pt() > pt_cut, "pt", p4.Pt());
    double eta_cut = period == analysis::Period::Run2017 ? cuts::hh_bbtautau_2017::TauTau::tauID::eta : cuts::H_tautau_2016::TauTau::tauID::eta;
    cut(std::abs(p4.Eta()) < eta_cut, "eta", p4.Eta());
    if(productionMode == ProductionMode::hh){
        const auto dmFinding = tau->tauID("decayModeFinding");
        cut(dmFinding > decayModeFinding, "oldDecayMode", dmFinding);
    }
    const auto packedLeadTauCand = dynamic_cast<const pat::PackedCandidate*>(tau->leadChargedHadrCand().get());
    cut(std::abs(packedLeadTauCand->dz()) < dz, "dz", packedLeadTauCand->dz());
    cut(std::abs(tau->charge()) == absCharge, "charge", tau->charge());
    if(productionMode == ProductionMode::hh) {
        cut(tau->tauID("againstElectronVLooseMVA6") > againstElectronVLooseMVA6, "againstElectron");
        cut(tau->tauID("againstMuonLoose3") > againstMuonLoose3, "againstMuon");
        if(period == analysis::Period::Run2017) {
            cut(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017") > 0.5, "VVLooseIso");
        }
    }
}

void TupleProducer_tauTau::FillEventTuple(const SelectionResults& selection)
{
    using Channel = analysis::Channel;
    using EventPart = ntuple::StorageMode::EventPart;

    BaseTupleProducer::FillEventTuple(selection, previous_selection.get());
    eventTuple().channelId = static_cast<int>(Channel::TauTau);

    ntuple::StorageMode storageMode(eventTuple().storageMode);
    eventTuple().storageMode = storageMode.Mode();

    BaseTupleProducer::FillTau(selection);

    eventTuple.Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TupleProducer_tauTau);
