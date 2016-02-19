/*! Definition of GenParticle class for analysis.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <set>

#include <TLorentzVector.h>

#include "TreeProduction/interface/GenParticle.h"

#include "Particles.h"
#include "AnalysisTools/Core/include/exception.h"

namespace analysis {

class GenParticle;

typedef std::vector< GenParticle > GenParticleVector;
typedef std::vector< const GenParticle* > GenParticlePtrVector;
typedef std::vector< GenParticlePtrVector > GenParticleVector2D;
typedef std::set< const GenParticle* > GenParticleSet;
typedef std::map<particles::ParticleCode, GenParticleSet> ParticleCodeMap;
typedef std::vector<size_t> GenParticleIndexVector;

class GenParticle {
public:
    size_t index;
    particles::PdgParticle pdg;
    particles::Status status;
    TLorentzVector momentum;
    TVector3 vertex;
    int charge;
    GenParticlePtrVector mothers;
    GenParticlePtrVector daughters;

public:
    GenParticle() {}
    GenParticle(size_t n, const ntuple::GenParticleVector& genParticles)
        : index(n)
    {
        if(n >= genParticles.size())
            throw std::runtime_error("GenParticles index is out of range.");
        const ntuple::GenParticle& ntupleGenParticle = genParticles.at(n);
//        if(ntupleGenParticle.Index != index)
//            throw std::runtime_error("bad index");
        pdg = particles::PdgParticle(ntupleGenParticle.PdgId);
        status = particles::NameProvider<particles::Status>::Convert(ntupleGenParticle.Status);
        momentum.SetPtEtaPhiM(ntupleGenParticle.pt, ntupleGenParticle.eta, ntupleGenParticle.phi, ntupleGenParticle.mass);
        vertex = TVector3(ntupleGenParticle.X, ntupleGenParticle.Y, ntupleGenParticle.Z);
        charge = ntupleGenParticle.Charge;
    }

    void Initialize(const ntuple::GenParticleVector& ntupleGenParticles, GenParticleVector& particles)
    {
        const ntuple::GenParticle& ntupleGenParticle = ntupleGenParticles.at(index);
        mothers.reserve(ntupleGenParticle.Mother_Indexes.size());
        //daughters.reserve(ntupleGenParticle.Daughter_Indexes.size());
        for (unsigned motherIndex : ntupleGenParticle.Mother_Indexes){
            GenParticle* mother = &particles.at(motherIndex);
            mothers.push_back(mother);
            mother->daughters.push_back(this);
        }
    }

    const GenParticle* GetHardInteractionOrigin() const
    {
        for(const GenParticle* origin = this; origin->pdg == pdg; origin = origin->mothers.front()) {
            if(origin->status == particles::HardInteractionProduct)
                return origin;
            if(origin->mothers.size() != 1) break;
        }
        throw exception("Hard interaction origin not found for particle ") << pdg.Name() << ", index = " << index;
    }
};

class VisibleGenObject {
public:
    const GenParticle* origin;

    GenParticleSet finalStateChargedLeptons;
    GenParticleSet finalStateChargedHadrons;
    GenParticleSet finalStateNeutralHadrons;

    TLorentzVector chargedLeptonsMomentum;
    TLorentzVector chargedHadronsMomentum;
    TLorentzVector neutralHadronsMomentum;
    TLorentzVector visibleMomentum;
    TLorentzVector invisibleMomentum;

    GenParticleSet particlesProcessed;

public:

    explicit VisibleGenObject() : origin(nullptr) {}

    explicit VisibleGenObject(const GenParticle *_origin) : origin(_origin)
    {
        const GenParticleSet particlesToIgnore;
        CollectInfo(origin, particlesProcessed, particlesToIgnore);
    }

    VisibleGenObject(const GenParticle *_origin, const GenParticleSet& particlesToIgnore) : origin(_origin)
    {
        CollectInfo(origin, particlesProcessed, particlesToIgnore);
    }

    bool operator < (const VisibleGenObject& other) const
    {
        return origin < other.origin;
    }


private:
    void CollectInfo(const GenParticle* particle, GenParticleSet& particlesProcessed,
                     const GenParticleSet& particlesToIgnore)
    {
        if(particle->status == particles::Status::FinalStateParticle && particle->daughters.size() != 0)
            throw exception("Invalid gen particle");

        if(particlesProcessed.count(particle) || particlesToIgnore.count(particle)) return;
        particlesProcessed.insert(particle);
        if(particle->status == particles::Status::FinalStateParticle) {
            if(particles::neutrinos.count(particle->pdg.Code)) {
                invisibleMomentum += particle->momentum;
                return;
            }

            visibleMomentum += particle->momentum;
            if(particle->pdg.Code == particles::e || particle->pdg.Code == particles::mu) {
                finalStateChargedLeptons.insert(particle);
                chargedLeptonsMomentum += particle->momentum;
            } else if(particle->charge) {
                finalStateChargedHadrons.insert(particle);
                chargedHadronsMomentum += particle->momentum;
            } else {
                finalStateNeutralHadrons.insert(particle);
                neutralHadronsMomentum += particle->momentum;
            }

        }

        for(const GenParticle* daughter : particle->daughters){
            CollectInfo(daughter, particlesProcessed, particlesToIgnore);
        }
    }
};

typedef std::vector<VisibleGenObject> VisibleGenObjectVector;

class GenEvent {
public:
    GenParticleVector genParticles;
    ParticleCodeMap particleCodeMap;
    ParticleCodeMap hardParticleCodeMap;
    GenParticleSet primaryParticles;

public:
    void Initialize(const ntuple::GenParticleVector& ntupleGenParticles)
    {
        genParticles.clear();
        particleCodeMap.clear();
        hardParticleCodeMap.clear();
        primaryParticles.clear();
        genParticles.reserve(ntupleGenParticles.size());

        for(size_t n = 0; n < ntupleGenParticles.size(); ++n){
            const GenParticle particle(n, ntupleGenParticles);
            genParticles.push_back(particle);
        }
        for(auto& particle : genParticles){
            particle.Initialize(ntupleGenParticles, genParticles);
            if (!particle.mothers.size())
                primaryParticles.insert(&particle);

            if (particle.status == particles::Decayed_or_fragmented ||
                    particle.status == particles::FinalStateParticle)
                particleCodeMap[particle.pdg.Code].insert(&particle);
            else if(particle.status == particles::HardInteractionProduct)
                hardParticleCodeMap[particle.pdg.Code].insert(&particle);
        }
    }

    GenParticleSet GetParticles(const particles::ParticleCodes& particleCodes, double pt = 0) const
    {
        GenParticleSet results;
        for (const particles::ParticleCode& code : particleCodes){
            const ParticleCodeMap::const_iterator code_iter = particleCodeMap.find(code);
            if (code_iter == particleCodeMap.end())
                continue;
            for (const GenParticle* particle : code_iter->second){
                if (particle->momentum.Pt() <= pt) continue;
                results.insert(particle);
            }
        }
        return results;
    }

    GenParticleSet GetHardParticles(const particles::ParticleCodes& particleCodes, double pt = 0) const
    {
        GenParticleSet results;
        for (const particles::ParticleCode& code : particleCodes){
            const ParticleCodeMap::const_iterator code_iter = hardParticleCodeMap.find(code);
            if (code_iter == hardParticleCodeMap.end())
                continue;
            for (const GenParticle* particle : code_iter->second){
                if (particle->momentum.Pt() <= pt) continue;
                results.insert(particle);
            }
        }
        return results;
    }

    void Print() const
    {
        for (const GenParticle* particle : primaryParticles) {
            PrintChain(particle);
        }
    }

    void PrintChain(const GenParticle* particle, unsigned iteration = 0) const
    {
        const particles::PdgParticle pdgParticle(particle->pdg);
        const particles::Status particleStatus = particles::NameProvider<particles::Status>::Convert(particle->status);
        const TLorentzVector genParticle_momentum = particle->momentum;
        for (unsigned n = 0; n < iteration; ++n)
            std::cout << "  ";
        std::cout << "index=" << particle->index << " name=" << pdgParticle << " status=" << particleStatus
                  <<  " pt= " << genParticle_momentum.Pt() <<"\n";
        for(unsigned n = 0; n < particle->daughters.size(); ++n) {
            const GenParticle* daughter = particle->daughters.at(n);
                PrintChain(daughter,iteration+1);
        }
    }
};

inline bool FindDecayProducts(const GenParticle& genParticle, const particles::ParticleCodes& particleCodes,
                              GenParticlePtrVector& decayProducts, bool ignoreFinalStateRadiation = false)
{
    if (genParticle.status != particles::Decayed_or_fragmented){
        throw std::runtime_error("particle type not supported");
    }

    decayProducts.clear();
    const GenParticlePtrVector* daughters = &genParticle.daughters;
    size_t expected_nDaughters = particleCodes.size();
    if (daughters->size() == 0){
        const GenParticle* originalParticle = &genParticle;
        while(originalParticle->status != particles::HardInteractionProduct) {
            if(originalParticle->mothers.size() != 1)
                throw std::runtime_error("more than one mother per a gen particle is not supported");
            const GenParticle* mother = originalParticle->mothers.front();
            if(mother->pdg != originalParticle->pdg
                    || (mother->status != particles::HardInteractionProduct && mother->daughters.size() != 1))
                throw std::runtime_error("particle is not part of the hard interaction final state radiation"
                                         " corrections chain");
            originalParticle = mother;
        }
        daughters = &originalParticle->daughters;
        ++expected_nDaughters;
    }
    size_t effective_nDaughters = daughters->size();
    if(ignoreFinalStateRadiation) {
        for(const GenParticle* daughter : *daughters) {
            if(daughter->pdg.Code == particles::gamma)
                --effective_nDaughters;
        }
    }
    if (effective_nDaughters != expected_nDaughters)
        return false;
    std::set<size_t> taken_daughters;
    for (const particles::ParticleCode& code : particleCodes){
        bool daughter_found = false;
        for (size_t n = 0; n < daughters->size(); ++n){
            if (taken_daughters.count(n))
                continue;
            if (code != daughters->at(n)->pdg.Code)
                continue;
            const GenParticle* daughter = daughters->at(n);
            if (daughter->status == particles::HardInteractionProduct){
                bool grandDaughter_found = false;
                for (const GenParticle* grandDaughter : daughter->daughters){
                    if (grandDaughter->pdg == daughter->pdg &&
                            (grandDaughter->status == particles::Decayed_or_fragmented ||
                             grandDaughter->status == particles::FinalStateParticle) ){
                        grandDaughter_found = true;
                        daughter = grandDaughter;
                        break;
                    }
                }
                if (!grandDaughter_found)
                    throw std::runtime_error("grand daughter not found");
            }

            decayProducts.push_back(daughter);
            taken_daughters.insert(n);
            daughter_found = true;
            break;
        }
        if (!daughter_found) return false;
    }
    return true;
}

inline bool FindDecayProducts2D(const GenParticlePtrVector& genParticles,
                                const particles::ParticleCodes2D& particleCodes2D, GenParticleVector2D& decayProducts2D,
                                GenParticleIndexVector& indexes, bool ignoreFinalStateRadiation = false)
{
    std::set<size_t> taken_genParticles;
    if (genParticles.size() != particleCodes2D.size())
        throw std::runtime_error("mismatched vector size of particles");
    decayProducts2D.clear();
    for(const particles::ParticleCodes& codes : particleCodes2D){
        bool particleFound = false;
        for (size_t n = 0; n < genParticles.size(); ++n ){
            if (taken_genParticles.count(n)) continue;
            const GenParticle& genParticle = *genParticles.at(n);
            GenParticlePtrVector decayProducts;
            if (!FindDecayProducts(genParticle,codes,decayProducts, ignoreFinalStateRadiation)) continue;
            particleFound = true;
            decayProducts2D.push_back(decayProducts);
            taken_genParticles.insert(n);
            indexes.push_back(n);
        }
        if (!particleFound) return false;
    }
    return true;
}

inline bool HasMatchWithMCParticle(const TLorentzVector& candidateMomentum, const GenParticle* genParticle, double deltaR)
{
    return genParticle && candidateMomentum.DeltaR(genParticle->momentum) < deltaR;
}

inline bool HasMatchWithMCObject(const TLorentzVector& candidateMomentum, const VisibleGenObject* genObject,
                                 double deltaR, bool useVisibleMomentum = false)
{
    if(!genObject) return false;
    const TLorentzVector& momentum = useVisibleMomentum ? genObject->visibleMomentum : genObject->origin->momentum;
    return candidateMomentum.DeltaR(momentum) < deltaR;
}

template<typename Container>
inline GenParticleSet FindMatchedParticles(const TLorentzVector& candidateMomentum,
                                           const Container& genParticles, double deltaR)
{
    GenParticleSet matchedGenParticles;
    for (const GenParticle* genParticle : genParticles){
        if (HasMatchWithMCParticle(candidateMomentum, genParticle, deltaR))
            matchedGenParticles.insert(genParticle);
    }
    return matchedGenParticles;
}

template<typename Container>
inline VisibleGenObjectVector FindMatchedObjects(const TLorentzVector& candidateMomentum,
                                                 const Container& genObjects, double deltaR)
{
    VisibleGenObjectVector matchedGenObjects;
    for (const VisibleGenObject& genObject : genObjects){
        if (HasMatchWithMCObject(candidateMomentum, &genObject, deltaR))
            matchedGenObjects.push_back(genObject);
    }
    return matchedGenObjects;
}

} // analysis
