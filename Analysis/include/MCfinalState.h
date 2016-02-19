/*! Definition of MCfinalState class for analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "GenParticle.h"

namespace analysis {

namespace finalState {

struct bbTauTau{

    const GenParticle* resonance;
    const GenParticle* Higgs_TauTau;
    const GenParticle* Higgs_BB;
    VisibleGenObjectVector b_jets;
    VisibleGenObjectVector taus;
    VisibleGenObjectVector hadronic_taus;

    bbTauTau() : resonance(nullptr), Higgs_TauTau(nullptr), Higgs_BB(nullptr){}

    virtual ~bbTauTau() {}
    virtual void Reset()
    {
        resonance = Higgs_TauTau = Higgs_BB = nullptr;
        b_jets.clear();
        taus.clear();
        hadronic_taus.clear();
    }
};

struct bbETaujet : public bbTauTau {
    const GenParticle* electron;
    const VisibleGenObject* tau_jet;

    bbETaujet() : electron(nullptr), tau_jet(nullptr) {}
    virtual void Reset() override
    {
        bbTauTau::Reset();
        electron = nullptr;
        tau_jet = nullptr;
    }
};

struct bbMuTaujet : public bbTauTau{
    const GenParticle* muon;
    const VisibleGenObject* tau_jet;

    bbMuTaujet() : muon(nullptr), tau_jet(nullptr) {}
    virtual void Reset() override
    {
        bbTauTau::Reset();
        muon = nullptr;
        tau_jet = nullptr;
    }
};

struct bbTaujetTaujet : public bbTauTau {
    const VisibleGenObject* leading_tau_jet;
    const VisibleGenObject* subleading_tau_jet;

    bbTaujetTaujet() : leading_tau_jet(nullptr), subleading_tau_jet(nullptr) {}
    virtual void Reset() override
    {
        bbTauTau::Reset();
        leading_tau_jet = subleading_tau_jet = nullptr;
    }
};

struct TauTau {
    const GenParticle* resonance;
    GenParticlePtrVector taus;
    GenParticlePtrVector hadronic_taus;
    TauTau() : resonance(nullptr) {}
    virtual ~TauTau() {}
};

struct ETaujet : public TauTau {
    const GenParticle* electron;
    const GenParticle* tau_jet;

    ETaujet() : electron(nullptr), tau_jet(nullptr) {}
    ETaujet(const TauTau& other) : TauTau(other), electron(nullptr), tau_jet(nullptr) {}
};

struct MuTaujet : public TauTau{
    const GenParticle* muon;
    const GenParticle* tau_jet;

    MuTaujet() : muon(nullptr), tau_jet(nullptr) {}
    MuTaujet(const TauTau& other) : TauTau(other), muon(nullptr), tau_jet(nullptr) {}
};

struct TaujetTaujet : public TauTau {
    const GenParticle* leading_tau_jet;
    const GenParticle* subleading_tau_jet;

    TaujetTaujet() : leading_tau_jet(nullptr), subleading_tau_jet(nullptr) {}
    TaujetTaujet(const TauTau& other) : TauTau(other), leading_tau_jet(nullptr), subleading_tau_jet(nullptr) {}
};

} // finalState
} // analysis
