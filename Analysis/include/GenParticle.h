/*! Definition of GenParticle class
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <set>
#include "h-tautau/Core/include/EventTuple.h"
#include "h-tautau/Analysis/include/GenStatusFlags.h"
#include "h-tautau/Analysis/include/Particle.h"
#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/exception.h"
#include <fstream>
namespace analysis {

class GenParticle;

using GenParticleVector = std::vector<GenParticle>;
using GenParticleSet = std::set<const GenParticle*>;
using GenParticlePtrVector = std::vector<const GenParticle*>;
using ParticleCodeMap = std::map<int,GenParticleSet>;

class GenParticle {
public:
    size_t index;
    int pdg;
    int status;
    GenStatusFlags genStatusFlags;
    LorentzVectorM momentum;
    GenParticlePtrVector mothers;
    GenParticlePtrVector daughters;
    Point3D vertex;

public:
    GenParticle(const ntuple::Event& events, size_t n);
};

class GenEvent {
private:
    GenParticleVector genParticles;
    ParticleCodeMap particleCodeMap;
    GenParticleSet primaryParticles;

public:
    GenEvent(const ntuple::Event& event);

    GenParticleSet GetParticles(int particle_pgd, bool requireIsLastCopy) const;

    void GetChosenParticlesTypes(const std::set<particles::ParticleType>& type_names, const GenParticle* mother, GenParticleSet& result) const;

    bool areParented(const GenParticle* daughter, const GenParticle* mother) const;

    static const std::string& GetParticleName(int pdgId);

    static void InitializeParticleDataTable(const std::string& fileName);

    void PrintChain(const GenParticle* particle, const std::string& pre = "") const;

    void Print() const;

    void FindFinalStateDaughters(const GenParticle& particle, std::set<const GenParticle*>& daughters,
                                 const std::set<int>& pdg_to_exclude);
    LorentzVectorM GetFinalStateMomentum(const GenParticle& particle, std::vector<const GenParticle*>& visible_daughters,
                                       bool excludeInvisible, bool excludeLightLeptons);

    static int GetParticleCharge(int pdg);
    static particles::ParticleType GetParticleType(int pdg);

private:
    static const std::unique_ptr<std::map<int, std::string>> particle_names;
    static const std::unique_ptr<std::map<int, particles::ParticleType>> particle_types;
    static const std::unique_ptr<std::map<int, int>> particle_charge;

};
} //analysis
