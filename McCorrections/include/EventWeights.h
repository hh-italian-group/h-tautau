/*! Various event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "PileUpWeight.h"
#include "LeptonWeights.h"
#include "BTagWeight.h"
#include "TopPtWeight.h"

namespace analysis {
namespace mc_corrections {

class EventWeights {
public:
    using Event = ntuple::Event;
    using ExpressEvent = ntuple::ExpressEvent;
    using PileUpWeightPtr = std::shared_ptr<PileUpWeight>;
    using LeptonWeightsPtr = std::shared_ptr<LeptonWeights>;
    using BTagWeightPtr = std::shared_ptr<BTagWeight>;
	using TopPtWeightPtr = std::shared_ptr<TopPtWeight>;

    EventWeights(Period period, DiscriminatorWP btag_wp)
    {
        if(period == Period::Run2015) {
            pileUp = PileUpWeightPtr(new class PileUpWeight(
                FullName("reWeight_Fall.root"), "lumiWeights", 60, 0));
            lepton = LeptonWeightsPtr(new LeptonWeights(
                FullLeptonName("Electron/Electron_IdIso0p10_eff.root"),
                FullLeptonName("Electron/Electron_SingleEle_eff.root"),
                FullLeptonName("Muon/Muon_IdIso0p1_fall15.root"),
                FullLeptonName("Muon/Muon_IsoMu18_fall15.root")));
            bTag = BTagWeightPtr(new class BTagWeight(
                FullName("bTagEff_Loose.root"),
                FullName("CSVv2.csv"),
                btag_wp));
        }
        else if(period == Period::Run2016) {
//            pileUp = PileUpWeightPtr(new class PileUpWeight(
//                FullName("Data_Pileup_2016_271036-276811_13TeVSpring16_PromptReco_69p2mbMinBiasXS.root"),
//                "pileup", 40, 0));
            pileUp = PileUpWeightPtr(new PileUpWeight(FullName("purewHisto.root"), "pileup", 60, 0));

            lepton = LeptonWeightsPtr(new LeptonWeights(
                FullLeptonName("Electron/Run2016BCD/Electron_IdIso0p10_eff.root"),
                FullLeptonName("Electron/Run2016BCD/Electron_Ele25eta2p1WPTight_eff.root"),
                FullLeptonName("Muon/Run2016BCD/Muon_IdIso0p15_eff.root"),
                FullLeptonName("Muon/Run2016BCD/Muon_IsoMu22_eff.root")));

            bTag = BTagWeightPtr(new class BTagWeight(
                FullName("bTagEfficiencies_80X.root"),
                FullName("CSVv2_ichep.csv"),
                btag_wp));
			
			top = TopPtWeightPtr(new TopPtWeight(0.0615, 0.0005));
			

        } else {
            throw exception("Period %1% is not supported.") % period;
        }
    }

	template<typename Event>
	double GetPileUpWeight(const Event& event) const { return pileUp ? pileUp->Get(event) : 1.; }
    double GetLeptonIdIsoWeight(const Event& event) const { return lepton ? lepton->GetIdIsoWeight(event) : 1.; }
    double GetLeptonTriggerWeight(const Event& event) const { return lepton ? lepton->GetTriggerWeight(event) : 1.; }
    double GetLeptonTotalWeight(const Event& event) const { return lepton ? lepton->GetTotalWeight(event) : 1.; }
    double GetBtagWeight(const Event& event) const { return bTag ? bTag->Compute(event) : 1.; }

	double GetTopPtWeight(const Event& event) const {return top ? top->Get(event) : 1.; }
    double GetTopPtWeight_express(const ExpressEvent& event) const {return top ? top->Get_express(event) : 1.; }

    double GetTotalWeight(const Event& event, bool apply_btag_weight = false)
    {
        double weight = GetPileUpWeight(event) * GetLeptonTotalWeight(event) * GetTopPtWeight(event);
        if(apply_btag_weight)
            weight *= GetBtagWeight(event);
        return weight;
    }

private:
    static std::string FullName(const std::string& fileName, const std::string& path)
    {
        return path + "/" + fileName;
    }

    static std::string FullName(const std::string& fileName)
    {
        static const std::string path = "h-tautau/McCorrections/data";
        return FullName(fileName, path);
    }

    static std::string FullLeptonName(const std::string& fileName)
    {
        static const std::string path = "HTT-utilities/LepEffInterface/data";
        return FullName(fileName, path);
    }

private:
    PileUpWeightPtr pileUp;
    LeptonWeightsPtr lepton;
    BTagWeightPtr bTag;
	TopPtWeightPtr top;
};

} // namespace mc_corrections
} // namespace analysis
