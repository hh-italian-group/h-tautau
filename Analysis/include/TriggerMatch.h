/*! Common tools and definitions suitable for analysis purposes.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <vector>
#include <string>
#include <memory>

#include <TH1.h>
#include <TLorentzVector.h>

#include "Particles.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "GenParticle.h"
#include "Candidate.h"

namespace analysis {

static const particles::ParticleCodes TauMuonicDecay = { particles::mu, particles::nu_mu, particles::nu_tau };
static const particles::ParticleCodes TauElectronDecay = { particles::e, particles::nu_e, particles::nu_tau };

inline bool IsLeptonicTau(const GenParticle& genParticle)
{
    if (genParticle.pdg.Code != particles::tau)
        throw exception("genParticle is not a tau!");
    GenParticlePtrVector tauProducts;
    return FindDecayProducts(genParticle,TauMuonicDecay,tauProducts,false) ||
                FindDecayProducts(genParticle,TauElectronDecay,tauProducts,false);
}

inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
                               const std::set<std::string>& interestinghltPaths, size_t& n)
{
    for (; n < objectMatchedPaths.size(); ++n){
        for (const std::string& interestingPath : interestinghltPaths){
            const std::string& objectMatchedPath = objectMatchedPaths.at(n);
            const size_t found = objectMatchedPath.find(interestingPath);
            if (found != std::string::npos) return true;
        }
    }
    return false;
}

inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
                               const std::set<std::string>& interestinghltPaths)
{
    size_t n = 0;
    return HaveTriggerMatched(objectMatchedPaths, interestinghltPaths, n);
}

inline bool HaveTriggerMatched(const ntuple::TriggerObjectVector& triggerObjects, const std::string& interestingPath,
                               const Candidate& candidate, double deltaR_Limit)
{
    if(candidate.GetFinalStateDaughters().size()) {
        for(const auto& daughter : candidate.GetFinalStateDaughters()) {
            if(!HaveTriggerMatched(triggerObjects, interestingPath, *daughter, deltaR_Limit))
                return false;
        }
        return true;
    }

    for (const ntuple::TriggerObject& triggerObject : triggerObjects){
        TLorentzVector triggerObjectMomentum;
        triggerObjectMomentum.SetPtEtaPhiM(triggerObject.pt, triggerObject.eta, triggerObject.phi, triggerObject.mass);
        for (unsigned n = 0; n < triggerObject.pathNames.size(); ++n){
            const std::string& objectMatchedPath = triggerObject.pathNames.at(n);
            const size_t found = objectMatchedPath.find(interestingPath);
            if (found != std::string::npos && triggerObject.pathValues.at(n) == 1 &&
                    triggerObjectMomentum.DeltaR(candidate.GetMomentum()) < deltaR_Limit)
                return true;
        }
    }
    return false;
}

inline bool HaveTriggerMatched(const std::string& interestingPath, const Candidate& candidate)
{
    if(candidate.GetFinalStateDaughters().size()) {
        for(const auto& daughter : candidate.GetFinalStateDaughters()) {
            if(!HaveTriggerMatched(interestingPath, *daughter))
                return false;
        }
        return true;
    }

    std::vector<std::string> objectMatchedPaths;
    if(candidate.GetType() == Candidate::Type::Tau)
        objectMatchedPaths = candidate.GetNtupleObject<ntuple::Tau>().matchedTriggerPaths;
    else if(candidate.GetType() == Candidate::Type::Jet)
        objectMatchedPaths = candidate.GetNtupleObject<ntuple::Jet>().matchedTriggerPaths;
    else if(candidate.GetType() == Candidate::Type::Muon)
        objectMatchedPaths = candidate.GetNtupleObject<ntuple::Muon>().matchedTriggerPaths;
    else if(candidate.GetType() == Candidate::Type::Electron)
        objectMatchedPaths = candidate.GetNtupleObject<ntuple::Electron>().matchedTriggerPaths;
    else
        throw exception("Unknow candidate to match trigger.");

    for (unsigned n = 0; n < objectMatchedPaths.size(); ++n){
        const std::string& objectMatchedPath = objectMatchedPaths.at(n);
        const size_t found = objectMatchedPath.find(interestingPath);
        if (found != std::string::npos) return true;
    }
    return false;
}

inline bool HaveTriggerMatched(const ntuple::TriggerObjectVector& triggerObjects,
                               const std::set<std::string>& interestingPaths, const Candidate& candidate,
                               double deltaR_Limit)
{
    for (const std::string& interestinPath : interestingPaths){
        if (HaveTriggerMatched(triggerObjects,interestinPath,candidate, deltaR_Limit)) return true;
    }
    return false;
}

inline bool HaveTriggerMatched(const std::set<std::string>& interestingPaths, const Candidate& candidate)
{
    for (const std::string& interestinPath : interestingPaths){
        if (HaveTriggerMatched(interestinPath, candidate)) return true;
    }
    return false;
}

} // analysis
