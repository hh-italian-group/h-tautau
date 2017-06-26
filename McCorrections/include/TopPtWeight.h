/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

class TopPtWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;
    using ExpressEvent = ntuple::ExpressEvent;

    TopPtWeight(double p1, double p2) : _p1(p1), _p2(p2) {}

    virtual double Get(const Event& event) const override
	{
		double topWeight = 1.;
		for (size_t n = 0; n < event.genParticles_pdg.size(); ++n)
		{
			if (std::abs(event.genParticles_pdg.at(n)) != 6) continue;
			const double pt = event.genParticles_p4.at(n).pt();
			topWeight *= std::sqrt(std::exp(_p1 - _p2 * pt));
		}
		return topWeight;
	}

    virtual double Get(const ExpressEvent& expr_event) const override
    {
        const double pt_top = expr_event.gen_top_pt;
        const double pt_topBar = expr_event.gen_topBar_pt;
        const double sf_1 = std::sqrt(std::exp(_p1 - _p2 * pt_top));
        const double sf_2 = std::sqrt(std::exp(_p1 - _p2 * pt_topBar));
        return sf_1 * sf_2;
    }

private:
    double _p1, _p2;
};

} // namespace mc_corrections
} // namespace analysis
