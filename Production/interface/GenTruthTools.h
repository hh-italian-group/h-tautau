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
    int n_chargedParticles;
    int n_neutralParticles;
};

void FindFinalStateDaughters(const reco::GenParticle& particle, std::set<const reco::GenParticle*>& daughters,
                             const std::set<int>& pdg_to_exclude);

LorentzVectorXYZ GetFinalStateMomentum(const reco::GenParticle& particle, std::vector<const reco::GenParticle*>& visible_daughters,
                                       bool excludeInvisible, bool excludeLightLeptons);

LeptonMatchResult LeptonGenMatch(const LorentzVectorM& p4,
    const reco::GenParticleCollection& genParticles);


float GetNumberOfPileUpInteractions(edm::Handle<std::vector<PileupSummaryInfo>>& pu_infos);

struct LheSummary {
    size_t n_partons = 0, n_b_partons = 0, n_c_partons = 0;
    double HT = 0., m_H = 0., m_hh = 0., cosTheta_hh = 0.;
    std::vector<int> index;
    std::vector<int> pdgId;
    std::vector<int> first_mother_index;
    std::vector<int> last_mother_index;
    std::vector<LorentzVectorM> p4;

};

LheSummary ExtractLheSummary(const LHEEventProduct& lheEventProduct);

} // namespace gen_truth
} // namespace analysis
