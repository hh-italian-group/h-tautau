/*! Definition of PatProduction classes for which ROOT dictionary should be made.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/PatProduction/interface/PatVertex.h"

namespace {
struct dictionary {
    pat::Vertex vertex;
    std::vector<pat::Vertex>::const_iterator v_p_v_ci;
    edm::Wrapper<std::vector<pat::Vertex> >  w_v_p_v;
    std::vector<float> v_f;
    std::vector<math::XYZVector> v_3d;
};
}
