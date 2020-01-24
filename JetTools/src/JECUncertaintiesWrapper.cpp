/*! Apply jet uncertainties to the event.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#include "h-tautau/JetTools/include/JECUncertaintiesWrapper.h"

namespace jec {

const std::set<UncertaintySource>& JECUncertaintiesWrapper::JetFullUncertainties()
{
    static const std::set<UncertaintySource> jetUncertainties = {
        UncertaintySource::JetFull_AbsoluteStat,
        UncertaintySource::JetFull_AbsoluteScale,
        UncertaintySource::JetFull_AbsoluteMPFBias,
        UncertaintySource::JetFull_AbsoluteFlavMap,
        UncertaintySource::JetFull_Fragmentation,
        UncertaintySource::JetFull_SinglePionECAL,
        UncertaintySource::JetFull_SinglePionHCAL,
        UncertaintySource::JetFull_FlavorQCD,
        UncertaintySource::JetFull_FlavorZJet,
        UncertaintySource::JetFull_FlavorPhotonJet,
        UncertaintySource::JetFull_FlavorPureGluon,
        UncertaintySource::JetFull_FlavorPureQuark,
        UncertaintySource::JetFull_FlavorPureCharm,
        UncertaintySource::JetFull_FlavorPureBottom,
        UncertaintySource::JetFull_TimePtEta,
        UncertaintySource::JetFull_RelativeJEREC1,
        UncertaintySource::JetFull_RelativeJEREC2,
        UncertaintySource::JetFull_RelativeJERHF,
        UncertaintySource::JetFull_RelativePtBB,
        UncertaintySource::JetFull_RelativePtEC1,
        UncertaintySource::JetFull_RelativePtEC2,
        UncertaintySource::JetFull_RelativePtHF,
        UncertaintySource::JetFull_RelativeBal,
        UncertaintySource::JetFull_RelativeFSR,
        UncertaintySource::JetFull_PileUpDataMC,
        UncertaintySource::JetFull_PileUpPtRef,
        UncertaintySource::JetFull_PileUpPtBB,
        UncertaintySource::JetFull_PileUpPtEC1,
        UncertaintySource::JetFull_PileUpPtEC2,
        UncertaintySource::JetFull_PileUpPtHF,
        UncertaintySource::JetFull_SubTotalPileUp,
        UncertaintySource::JetFull_SubTotalRelative,
        UncertaintySource::JetFull_SubTotalPt,
        UncertaintySource::JetFull_SubTotalScale,
        UncertaintySource::JetFull_SubTotalAbsolute,
        UncertaintySource::JetFull_SubTotalMC,
        UncertaintySource::JetFull_TotalNoFlavor,
        UncertaintySource::JetFull_TotalNoTime,
        UncertaintySource::JetFull_TotalNoFlavorNoTime
    };
    return jetUncertainties;
}

const std::set<UncertaintySource>& JECUncertaintiesWrapper::JetReducedUncertainties()
{
    static const std::set<UncertaintySource> jetUncertainties = {
        UncertaintySource::JetReduced_Absolute,
        UncertaintySource::JetReduced_Absolute_year,
        UncertaintySource::JetReduced_BBEC1,
        UncertaintySource::JetReduced_BBEC1_year,
        UncertaintySource::JetReduced_EC2,
        UncertaintySource::JetReduced_EC2_year,
        UncertaintySource::JetReduced_FlavorQCD,
        UncertaintySource::JetReduced_HF,
        UncertaintySource::JetReduced_HF_year,
        UncertaintySource::JetReduced_RelativeBal,
        UncertaintySource::JetReduced_RelativeSample_year
    };
    return jetUncertainties;
}

bool JECUncertaintiesWrapper::IsJetUncertainties(UncertaintySource unc_source)
{
    if(JECUncertaintiesWrapper::JetReducedUncertainties().count(unc_source) ||
       JECUncertaintiesWrapper::JetFullUncertainties().count(unc_source) ||
       unc_source == UncertaintySource::JetFull_Total ||
       unc_source == UncertaintySource::JetReduced_Total)
        return true;
    return false;
}

const std::string JECUncertaintiesWrapper::ReturnJecName(UncertaintySource unc_source, bool is_full, analysis::Period& period)
{
    std::string new_string;
    if(is_full){
        std::string full_name = analysis::ToString(unc_source);
        std::string jet_string = "JetFull_";
        new_string = full_name.erase(0,jet_string.size() - 1);
    }
    else{
        std::string full_name = analysis::ToString(unc_source);
        std::string jet_string = "JetReduced_";
        std::string partial_string = full_name.erase(0,jet_string.size());
        std::string year = "year";
        size_t found = partial_string.find(year);
        if(found != std::string::npos){
            std::string period_name = analysis::ToString(period);
            std::string run_string = "Run";
            std::string period_string = period_name.erase(0,run_string.size());
            std::cout << "period_string: " << period_string <<std::endl;
            new_string = partial_string.replace(found,year.length(),period_string);
        }
        else{
            new_string = partial_string;
        }
    }
    std::cout << "new_string: " << new_string << std::endl;
    return new_string;
}

JECUncertaintiesWrapper::JECUncertaintiesWrapper(const std::string& uncertainties_source, bool is_full, analysis::Period& period)
{
    auto createUncSet = [](bool is_full) {
        std::set<UncertaintySource> jetUncertainties;
        if(is_full){
            jetUncertainties = JECUncertaintiesWrapper::JetFullUncertainties();
            jetUncertainties.insert(UncertaintySource::JetFull_Total);
        }
        else{
            jetUncertainties = JECUncertaintiesWrapper::JetReducedUncertainties();
            jetUncertainties.insert(UncertaintySource::JetReduced_Total);
        }
        return jetUncertainties;
     };

    static const std::set<UncertaintySource> jetUncertaintiesTotal = createUncSet(is_full);

    for (const auto jet_unc : jetUncertaintiesTotal) {
        std::string full_name = JECUncertaintiesWrapper::ReturnJecName(jet_unc,is_full,period);
        std::cout << "Full name: " << full_name << std::endl;
        JetCorrectorParameters p(uncertainties_source, full_name);
        auto unc = std::make_shared<JetCorrectionUncertainty>(p);
        uncertainty_map[jet_unc] = unc;
    }
}

}
