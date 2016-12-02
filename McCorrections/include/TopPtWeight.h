/*! The pile up weight.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"


namespace analysis {
namespace mc_corrections {

class TopPtWeight {

public:
	
	using Event = ntuple::Event;

	TopPtWeight(double p1, double p2):
	_p1(p1), _p2(p2) {}

	template<typename Event>
	double Get(const Event& event) const
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

private:
	double _p1;
	double _p2;

};

} // namespace mc_corrections
} // namespace analysis
