/*! Tools for working with MC generator truth.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"

namespace analysis {

namespace gen_truth {

using MatchResult = std::pair<GenMatch, const reco::GenParticle*>;

inline void FindFinalStateDaughters(const reco::GenParticle& particle, std::vector<const reco::GenParticle*>& daughters,
                                    const std::set<int>& pdg_to_exclude = {})
{
    if(!particle.daughterRefVector().size()) {
        const int abs_pdg = std::abs(particle.pdgId());
        if(!pdg_to_exclude.count(abs_pdg))
            daughters.push_back(&particle);
    } else {
        for(const auto& daughter : particle.daughterRefVector())
            FindFinalStateDaughters(*daughter, daughters, pdg_to_exclude);
    }
}

inline LorentzVectorXYZ GetFinalStateMomentum(const reco::GenParticle& particle, bool excludeInvisible,
                                              bool excludeLightLeptons)
{
    using set = std::set<int>;
    using pair = std::pair<bool, bool>;
    static const set empty = {};
    static const set light_leptons = { 11, 13 };
    static const set invisible_particles = { 12, 14, 16 };
    static const set light_and_invisible = tools::union_sets({light_leptons, invisible_particles});

    static const std::map<pair, const set*> to_exclude {
        { pair(false, false), &empty }, { pair(true, false), &invisible_particles },
        { pair(false, true), &light_leptons }, { pair(true, true), &light_and_invisible },
    };

    std::vector<const reco::GenParticle*> daughters;
    FindFinalStateDaughters(particle, daughters, *to_exclude.at(pair(excludeInvisible, false)));

    LorentzVectorXYZ p4;
    for(auto daughter : daughters){
	if(excludeLightLeptons && light_leptons.count(std::abs(daughter->pdgId())) && daughter->statusFlags().isDirectTauDecayProduct()) continue;     
	p4 += daughter->p4();
    }
    return p4;
}

template<typename LVector>
MatchResult LeptonGenMatch(const LVector& p4, const std::vector<reco::GenParticle>& genParticles)
{
    static constexpr int electronPdgId = 11, muonPdgId = 13, tauPdgId = 15;
    static constexpr double dR2_threshold = std::pow(0.2, 2);

    static const std::map<int, double> pt_thresholds = {
        { electronPdgId, 8 }, { muonPdgId, 8 }, { tauPdgId, 15 }
    };

    using pair = std::pair<int, bool>;
    static const std::map<pair, GenMatch> genMatches = {
        { { electronPdgId, false }, GenMatch::Electron }, { { electronPdgId, true }, GenMatch::TauElectron },
        { { muonPdgId, false }, GenMatch::Muon }, { { muonPdgId, true }, GenMatch::TauMuon },
        { { tauPdgId, false }, GenMatch::Tau }, { { tauPdgId, true }, GenMatch::Tau }
    };

    MatchResult result(GenMatch::NoMatch, nullptr);
    double match_dr2 = dR2_threshold;
 


    for(const reco::GenParticle& particle : genParticles) {
        const bool isTauProduct = particle.statusFlags().isDirectPromptTauDecayProduct();
        if((!particle.statusFlags().isPrompt() && !isTauProduct) || !particle.statusFlags().isLastCopy()) continue;

        const int abs_pdg = std::abs(particle.pdgId());
        if(!pt_thresholds.count(abs_pdg)) continue;

        const auto particle_p4 = abs_pdg == tauPdgId ? GetFinalStateMomentum(particle, true, true) : particle.p4();

        const double dr2 = ROOT::Math::VectorUtil::DeltaR2(p4, particle_p4);
        if(dr2 >= match_dr2) continue;
        if(particle_p4.pt() <= pt_thresholds.at(abs_pdg)) continue;

        match_dr2 = dr2;
        result.first = genMatches.at(pair(abs_pdg, isTauProduct));
        result.second = &particle;
    }
    return result;
}

} // namespace gen_truth
} // namespace analysis
