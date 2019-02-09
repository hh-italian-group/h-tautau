/*! Apply jet uncertainties to the event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"

namespace jec {

const std::set<UncertaintySource>& JECUncertaintiesWrapper::JetUncertainties()
{
    static const std::set<UncertaintySource> jetUncertainties = {
        UncertaintySource::AbsoluteStat,
        UncertaintySource::AbsoluteScale,
        UncertaintySource::AbsoluteMPFBias,
        UncertaintySource::AbsoluteFlavMap,
        UncertaintySource::Fragmentation,
        UncertaintySource::SinglePionECAL,
        UncertaintySource::SinglePionHCAL,
        UncertaintySource::FlavorQCD,
        UncertaintySource::FlavorZJet,
        UncertaintySource::FlavorPhotonJet,
        UncertaintySource::FlavorPureGluon,
        UncertaintySource::FlavorPureQuark,
        UncertaintySource::FlavorPureCharm,
        UncertaintySource::FlavorPureBottom,
        UncertaintySource::TimePtEta,
        UncertaintySource::RelativeJEREC1,
        UncertaintySource::RelativeJEREC2,
        UncertaintySource::RelativeJERHF,
        UncertaintySource::RelativePtBB,
        UncertaintySource::RelativePtEC1,
        UncertaintySource::RelativePtEC2,
        UncertaintySource::RelativePtHF,
        UncertaintySource::RelativeBal,
        UncertaintySource::RelativeFSR,
        UncertaintySource::PileUpDataMC,
        UncertaintySource::PileUpPtRef,
        UncertaintySource::PileUpPtBB,
        UncertaintySource::PileUpPtEC1,
        UncertaintySource::PileUpPtEC2,
        UncertaintySource::PileUpPtHF,
        UncertaintySource::SubTotalPileUp,
        UncertaintySource::SubTotalRelative,
        UncertaintySource::SubTotalPt,
        UncertaintySource::SubTotalScale,
        UncertaintySource::SubTotalAbsolute,
        UncertaintySource::SubTotalMC,
        UncertaintySource::TotalNoFlavor,
        UncertaintySource::TotalNoTime,
        UncertaintySource::TotalNoFlavorNoTime
    };
    return jetUncertainties;
}

const std::set<UncertaintySource>& JECUncertaintiesWrapper::JetUncertainties_withTotal()
{
    auto createUncSet = []() {
        std::set<UncertaintySource> jetUncertainties = JetUncertainties();
        jetUncertainties.insert(UncertaintySource::Total);
        return jetUncertainties;
     };

    static const std::set<UncertaintySource> jetUncertaintiesTotal = createUncSet();
    return jetUncertaintiesTotal;
}

JECUncertaintiesWrapper::JECUncertaintiesWrapper(const std::string& uncertainties_source)
{
    for (const auto jet_unc : JetUncertainties_withTotal()) {
        std::string full_name = analysis::ToString(jet_unc);
        JetCorrectorParameters p(uncertainties_source, full_name);
        auto unc = std::make_shared<JetCorrectionUncertainty>(p);
        uncertainty_map[jet_unc] = unc;
    }
}

}
