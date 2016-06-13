/*! Tools for working with MC generator truth.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

namespace analysis {

namespace gen_truth {

inline int genMatch(const reco::Candidate::LorentzVector& p4, const std::vector<reco::GenParticle>& genParticles)
{
  int match = 6;
  float dr = 0.2;
  float dr_new = 0.2;

  for ( std::vector<reco::GenParticle>::const_iterator genp = genParticles.begin();
    genp != genParticles.end(); ++genp )
    {
    dr_new = ROOT::Math::VectorUtil::DeltaR( p4, genp->p4());

    if (dr_new < dr){
      if (TMath::Abs(genp->pdgId()) == 11) {
          if (genp->pt() > 8) {
             if (genp->statusFlags().isPrompt()) {
                match = 1;
                dr = dr_new;
             }
             else if (genp->statusFlags().isDirectPromptTauDecayProduct()) {
                match = 3;
                dr = dr_new;
             }
          }
      }
      else if (TMath::Abs(genp->pdgId()) == 13) {
          if (genp->pt() > 8) {
             if (genp->isPromptFinalState()){
                match = 2;
                dr = dr_new;
             }
             if (genp->isDirectPromptTauDecayProductFinalState()){
                match = 4;
                dr = dr_new;
             }
          }
      }
      else if (TMath::Abs(genp->pdgId()) == 15) {
          if (genp->statusFlags().isPrompt()) {
             reco::Candidate::LorentzVector tau_p4 = utils_genMatch::getVisMomentumNoLep(&(*genp));
            if (tau_p4.pt() > 15) {
                 match = 5;
                dr = dr_new;
             }
          }
      }
    }
  } //GenParticle Loop

  return match;
}

} // namespace gen_truth
} // namespace analysis
