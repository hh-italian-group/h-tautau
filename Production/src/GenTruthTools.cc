/*! Tools for working with MC generator truth.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/Production/interface/GenTruthTools.h"

namespace analysis {

namespace gen_truth {

void FindFinalStateDaughters(const reco::GenParticle& particle, std::set<const reco::GenParticle*>& daughters,
                             const std::set<int>& pdg_to_exclude)
{
    if(!particle.daughterRefVector().size()) {
        const int abs_pdg = std::abs(particle.pdgId());
        if(!pdg_to_exclude.count(abs_pdg))
            daughters.insert(&particle);
    } else {
        for(const auto& daughter : particle.daughterRefVector())
            FindFinalStateDaughters(*daughter, daughters, pdg_to_exclude);
    }
}

LorentzVectorXYZ GetFinalStateMomentum(const reco::GenParticle& particle, std::vector<const reco::GenParticle*>& visible_daughters,
                                       bool excludeInvisible, bool excludeLightLeptons)
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



    std::set<const reco::GenParticle*> daughters_set;
    FindFinalStateDaughters(particle, daughters_set, *to_exclude.at(pair(excludeInvisible, false)));
    visible_daughters.clear();
    visible_daughters.insert(visible_daughters.begin(), daughters_set.begin(), daughters_set.end());

    LorentzVectorXYZ p4;
    for(auto daughter : visible_daughters) {
        if(excludeLightLeptons && light_leptons.count(std::abs(daughter->pdgId()))
                && daughter->statusFlags().isDirectTauDecayProduct()) continue;
            p4 += daughter->p4();
    }
    return p4;
}

LeptonMatchResult LeptonGenMatch(const LorentzVectorM& p4, const reco::GenParticleCollection& genParticles)
{
    static constexpr int electronPdgId = 11, muonPdgId = 13, tauPdgId = 15;
    static double dR2_threshold = std::pow(0.2, 2);

    static const std::map<int, double> pt_thresholds = {
        { electronPdgId, 8 }, { muonPdgId, 8 }, { tauPdgId, 15 }
    };

    using pair = std::pair<int, bool>;
    static const std::map<pair, GenLeptonMatch> genMatches = {
        { { electronPdgId, false }, GenLeptonMatch::Electron }, { { electronPdgId, true }, GenLeptonMatch::TauElectron },
        { { muonPdgId, false }, GenLeptonMatch::Muon }, { { muonPdgId, true }, GenLeptonMatch::TauMuon },
        { { tauPdgId, false }, GenLeptonMatch::Tau }, { { tauPdgId, true }, GenLeptonMatch::Tau }
    };

    LeptonMatchResult result;
    double match_dr2 = dR2_threshold;

    for(const reco::GenParticle& particle : genParticles) {
        const bool isTauProduct = particle.statusFlags().isDirectPromptTauDecayProduct();
        if((!particle.statusFlags().isPrompt() && !isTauProduct) /*|| !particle.statusFlags().isLastCopy()*/) continue;

        const int abs_pdg = std::abs(particle.pdgId());
        if(!pt_thresholds.count(abs_pdg)) continue;

        std::vector<const reco::GenParticle*> visible_daughters;
        const auto particle_p4 = abs_pdg == tauPdgId ? GetFinalStateMomentum(particle, visible_daughters, true, true)
                                                     : particle.p4();

        const double dr2 = ROOT::Math::VectorUtil::DeltaR2(p4, particle_p4);
        if(dr2 >= match_dr2) continue;
        if(particle_p4.pt() <= pt_thresholds.at(abs_pdg)) continue;

        match_dr2 = dr2;
        result.match = genMatches.at(pair(abs_pdg, isTauProduct));
        result.gen_particle = &particle;
        result.visible_daughters = visible_daughters;
        result.visible_daughters_p4 = particle_p4;

        int n_chargedParticles = 0;
        int n_neutralParticles = 0;
        for(unsigned n = 0; n < visible_daughters.size(); ++n){
          const reco::GenParticle* gen_visible_particle = visible_daughters.at(n);
          if(gen_visible_particle->charge() == 0)
            ++n_neutralParticles;
          if(gen_visible_particle->charge() != 0)
            ++n_chargedParticles;
        }

        result.n_chargedParticles = n_chargedParticles;
        result.n_neutralParticles = n_neutralParticles;

    }
    return result;
}

} // namespace gen_truth
} // namespace analysis
