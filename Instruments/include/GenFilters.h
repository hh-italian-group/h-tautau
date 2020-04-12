/*! Definition of a various gen filters.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <vector>
#include <set>
#include <map>
#include <Math/PtEtaPhiE4D.h>
#include <Math/PtEtaPhiM4D.h>
#include <Math/LorentzVector.h>
#include <ROOT/RVec.hxx>

#pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > >;
#pragma link C++ class vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > >;

inline bool FilterHTT(const ROOT::VecOps::RVec<Int_t>& genParticles_index,
                      const ROOT::VecOps::RVec<Int_t>& genParticles_pdg,
                      const ROOT::VecOps::RVec<Int_t>& genParticles_rel_pIndex,
                      const ROOT::VecOps::RVec<Int_t>& genParticles_rel_mIndex)
{
    static constexpr int tauPdgId = 15, higgsPdgId = 25;
    std::set<int> taus;
    std::map<int, size_t> higgses;
    for(size_t n = 0; n < genParticles_pdg.size(); ++n) {
        const int abs_pdg = std::abs(genParticles_pdg.at(n));
        if(abs_pdg == higgsPdgId)
            higgses[genParticles_index.at(n)] = 0;
        else if(abs_pdg == tauPdgId)
            taus.insert(genParticles_index.at(n));
    }
    for(size_t n = 0; n < genParticles_rel_pIndex.size(); ++n) {
        if(taus.count(genParticles_rel_pIndex.at(n))) {
            auto iter = higgses.find(genParticles_rel_mIndex.at(n));
            if(iter != higgses.end()) {
                if(iter->second > 0)
                    return true;
                ++iter->second;
            }
        }
    }

    return false;
}
