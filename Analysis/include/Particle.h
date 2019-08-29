/*! Definiton of Particle
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <string>
#include <set>
#include "AnalysisTools/Core/include/Tools.h"
#include "AnalysisTools/Core/include/EnumNameMap.h"

namespace particles {

using ::analysis::operator<<;
using ::analysis::operator>>;

namespace ParticleCode {

    enum pdg {e= 11, mu =13, tau = 15 , b = 5, higgs = 25, nu_e = 12, nu_mu = 14, nu_tau = 16, gamma= 22, g= 21};

    ENUM_NAMES(pdg) = {
        { pdg::e, "e" },
        { pdg::mu, "mu" },
        { pdg::tau, "tau" },
        { pdg::b, "b" },
        { pdg::higgs, "higgs" },
        { pdg::nu_e, "nu_e" },
        { pdg::nu_mu, "nu_mu" },
        { pdg::nu_tau, "nu_tau" },
        { pdg::gamma, "gamma" },
        { pdg::g, "g" }
    };
};

enum class ParticleType {baryon, meson, lepton, diquark, boson, quark, nucleus};

ENUM_NAMES(ParticleType) = {
    { ParticleType::baryon, "baryon" },
    { ParticleType::meson, "meson" },
    { ParticleType::lepton, "lepton" },
    { ParticleType::diquark, "diquark"},
    { ParticleType::boson, "boson"},
    { ParticleType::quark, "quark"},
    { ParticleType::nucleus, "nucleus"}
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
