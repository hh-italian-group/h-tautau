/*! Tools for working with MC generator truth.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"


namespace analysis {

namespace gen_truth {



struct LeptonMatchResult {
    GenLeptonMatch match{GenLeptonMatch::NoMatch};
    const reco::GenParticle* gen_particle{nullptr};
    std::vector<const reco::GenParticle*> visible_daughters;
    LorentzVectorXYZ visible_daughters_p4;
};

void FindFinalStateDaughters(const reco::GenParticle& particle, std::set<const reco::GenParticle*>& daughters,
                             const std::set<int>& pdg_to_exclude);

LorentzVectorXYZ GetFinalStateMomentum(const reco::GenParticle& particle, std::vector<const reco::GenParticle*>& visible_daughters,
                                       bool excludeInvisible, bool excludeLightLeptons);

LeptonMatchResult LeptonGenMatch(const LorentzVectorM& p4,
    const reco::GenParticleCollection& genParticles);


inline float GetNumberOfPileUpInteractions(edm::Handle<std::vector<PileupSummaryInfo>>& pu_infos)
{
    if(pu_infos.isValid()) {
        for(const PileupSummaryInfo& pu : *pu_infos) {
            if(pu.getBunchCrossing() == 0)
                return pu.getTrueNumInteractions();
        }
    }
    return std::numeric_limits<float>::lowest();
}

struct LheSummary {
    size_t n_partons = 0, n_b_partons = 0, n_c_partons = 0;
    double HT = 0., m_H = 0., m_hh = 0., cosTheta_hh = 0.;
};

inline LheSummary ExtractLheSummary(const LHEEventProduct& lheEventProduct)
{
    static constexpr int c_quark = 4, b_quark = 5, H0 = 35, h0 = 25;
    static const std::set<int> quarks_and_gluons = { 1, 2, 3, 4, 5, 6, 21 };

    LheSummary summary;
    const lhef::HEPEUP& lheEvent = lheEventProduct.hepeup();
    const std::vector<lhef::HEPEUP::FiveVector>& lheParticles = lheEvent.PUP;
    std::vector<analysis::LorentzVectorXYZ> h0_p4;
    for(size_t n = 0; n < lheParticles.size(); ++n) {
        const int absPdgId = std::abs(lheEvent.IDUP[n]);
        const int status = lheEvent.ISTUP[n];
        if(absPdgId == H0) summary.m_H = lheParticles[n][4];
        if(absPdgId == h0)
            h0_p4.push_back(analysis::LorentzVectorXYZ(lheParticles[n][0], lheParticles[n][1],
                                                              lheParticles[n][2], lheParticles[n][3]));
        if(status != 1 || !quarks_and_gluons.count(absPdgId)) continue;
        ++summary.n_partons;
        if(absPdgId == c_quark) ++summary.n_c_partons;
        if(absPdgId == b_quark) ++summary.n_b_partons;
        summary.HT += std::sqrt(std::pow(lheParticles[n][0], 2) + std::pow(lheParticles[n][1], 2));
    }
    if(h0_p4.size() == 2) {
        const analysis::LorentzVectorXYZ H_p4 = h0_p4.at(0) + h0_p4.at(1);
        const auto boosted_h0 = ROOT::Math::VectorUtil::boost(h0_p4.at(0), H_p4.BoostToCM());
        summary.cosTheta_hh = ROOT::Math::VectorUtil::CosTheta(boosted_h0, ROOT::Math::Cartesian3D<>(0, 0, 1));
        summary.m_hh = H_p4.mass();
    }
    return summary;
}

} // namespace gen_truth
} // namespace analysis
