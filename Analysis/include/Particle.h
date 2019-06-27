/*! Definiton of Particle
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <string>
#include <set>
#include "AnalysisTools/Core/include/Tools.h"

namespace particles {
namespace ParticleCode {

    enum pdg {e= 11, mu =13, tau = 15 , b = 5, higgs = 25, nu_e = 12, nu_mu = 14, nu_tau = 16, gamma= 22, g= 21};
    enum ParticleType {baryon, meson, lepton, diquark};
};

inline const std::set<int>& neutrinos()
{
    static const std::set<int> neutrinos = { ParticleCode::nu_e, ParticleCode::nu_mu, ParticleCode::nu_tau };
    return neutrinos;
}

inline const std::set<int>& light_leptons()
{
    static const std::set<int> light_leptons = { ParticleCode::e, ParticleCode::mu };
    return light_leptons;
}

inline const std::set<int>& light_and_invisible()
{
//    auto neutrinos_ = neutrinos();
    const std::set<int>& neutrinos_ = neutrinos();
    const std::set<int>& light_leptons_ = light_leptons();
    static const std::set<int> light_and_invisible = analysis::tools::union_sets({light_leptons_, neutrinos_});
    return light_and_invisible;
}

} // particles
