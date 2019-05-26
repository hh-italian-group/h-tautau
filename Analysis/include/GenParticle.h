/*! Definition of data and event categories used in HH->bbTauTau analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <set>
#include <TLorentzVector.h>

#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Analysis/include/GenStatusFlags.h"
#include "h-tautau/Analysis/include/Particle.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/exception.h"
#include <iostream>
#include <fstream>

namespace analysis {

class GenParticle;

using GenParticleVector = std::vector<GenParticle>;
using GenParticleSet = std::set<const GenParticle*>;
using GenParticlePtrVector = std::vector<const GenParticle*>;
using ParticleCodeMap = std::map<int,GenParticleSet>;
using GenParticleVector2D = std::vector<GenParticlePtrVector>;


class GenParticle {
public:
    size_t index;
    int pdg;
    int status;
    GenStatusFlags genStatusFlags;
    LorentzVectorM_Float momentum;
    GenParticlePtrVector mothers;
    GenParticlePtrVector daughters;

public:
    GenParticle(const ntuple::Event& events, size_t n);
};

class GenEvent {
public:
    GenParticleVector genParticles;
    ParticleCodeMap particleCodeMap;
    GenParticleSet primaryParticles;

public:
    GenEvent(const ntuple::Event& event);

    GenParticleSet GetParticles(int particle_pgd) const;

    std::vector<const GenParticle*> GetTypesParticles(std::vector<std::string> type_names, const GenParticle* possible_mother ) const;

    bool areParented(const GenParticle* daughter, const GenParticle* possible_mother) const;

    static std::string& GetParticleName(int pdgId);

    static void intializeNames(const std::string& fileName);

    void PrintChain(const GenParticle* particle, const std::string& pre) const;

    void Print() const;

    void FindFinalStateDaughters(const GenParticle& particle, std::set<const GenParticle*>& daughters,
                                 const std::set<int>& pdg_to_exclude);
    LorentzVectorM_Float GetFinalStateMomentum(const GenParticle& particle, std::vector<const GenParticle*>& visible_daughters,
                                       bool excludeInvisible, bool excludeLightLeptons);



private:
    static std::map<int, std::string> particle_names;
    static std::map<int, std::string> particle_types;
};
} //analysis
