/*! Definition of wrappers for KinFit.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Core/include/SummaryTuple.h"
#include "h-tautau/Core/include/TupleObjects.h"
#include "h-tautau/Cuts/include/hh_bbtautau_2017.h"
#include "h-tautau/Cuts/include/H_tautau_2016_baseline.h"
#include "h-tautau/Cuts/include/H_tautau_2017_baseline.h"

namespace analysis {

enum class SignalMode { HTT, HTT_sync, TauPOG_default, TauPOG_deepTauVsJet, TauPOG_deepTauVsJet_full, TauPOG_dpfTau, HH, bbtautau, Skimmer, TauPOG_Skimmer };

ENUM_NAMES(SignalMode) = {
    { SignalMode::HTT, "HTT" },
    { SignalMode::HTT_sync, "HTT_sync" },
    { SignalMode::TauPOG_default, "TauPOG_default" },
    { SignalMode::TauPOG_deepTauVsJet, "TauPOG_deepTauVsJet" },
    { SignalMode::TauPOG_deepTauVsJet_full, "TauPOG_deepTauVsJet_full" },
    { SignalMode::HH, "HH" },
    { SignalMode::bbtautau, "bbtautau" },
    { SignalMode::Skimmer, "Skimmer" },
    { SignalMode::TauPOG_Skimmer, "TauPOG_Skimmer" }
};

class SignalObjectSelector {
public:
    SignalObjectSelector(SignalMode _mode);

    bool PassLeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const;
    boost::optional<size_t> GetHiggsCandidateIndex(const ntuple::Event& event) const;
    bool PassLeptonVetoSelection(const ntuple::Event& event) const;
    bool PassMETfilters(const ntuple::Event& event) const;


private:
    bool PassHTT_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel, bool is_sync) const;
    bool PassTauPOG_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const;
    bool PassHH_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const;
    bool Pass_bbtautau_LeptonSelection(const ntuple::TupleLepton& lepton, Channel channel) const;
    bool PassSkimmer_LeptonSelection(const ntuple::TupleLepton& lepton) const;
    bool PassTauPOG_Skimmer_LeptonSelection(const ntuple::TupleLepton& lepton) const;


private:
    SignalMode mode;
    double DR2_leptons;
    TauIdDiscriminator discriminator;
};

} // namespace analysis
