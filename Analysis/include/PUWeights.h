/*!
 * \file EventWeights.h
 * \brief Definition of the class to calculate and store different event weights.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-11-17 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "Tools.h"
#include "SyncTree.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"

namespace analysis {

class PUWeights {
public:
    PUWeights(bool _is_data, bool _is_embedded, bool _apply_pu_weight, bool _apply_DM_weight,
                 const std::string& pu_reweight_file_name, double _max_available_pu, double _default_pu_weight)
        : is_data(_is_data), is_embedded(_is_embedded), apply_pu_weight(_apply_pu_weight),
          max_available_pu(_max_available_pu),
          default_pu_weight(_default_pu_weight)
    {

        isoIdScaleFactor = new htt_utilities::ScaleFactor();
        triggerScaleFactor = new htt_utilities::ScaleFactor();
        isoIdScaleFactor->init_ScaleFactor("/Users/claup/workspace/HHbbTauTau_Run2/HTT-utilities/LepEffInterface/data/Muon/Muon_IdIso0p1_fall15.root");
        triggerScaleFactor->init_ScaleFactor("/Users/claup/workspace/HHbbTauTau_Run2/HTT-utilities/LepEffInterface/data/Muon/Muon_IsoMu18_fall15.root");

        if(is_data && apply_pu_weight)
            throw exception("Inconsistend event weight configuration: requested to apply PU weight for data sample.");
//        if(is_embedded && apply_pu_weight)
//            throw exception("Inconsistend event weight configuration: requested to apply PU weight for embedded sample.");
        if(is_data && is_embedded)
            throw exception("Inconsistend event weight configuration: sample is data and embedded at the same time.");
        if(apply_pu_weight)
            pu_weights = LoadPUWeights(pu_reweight_file_name);
        Reset();
    }

    virtual void Reset()
    {
        eventWeight = 1;
        PUweight = 1;
        isoIdLeptonWeight = 1;
        triggerLeptonWeight = 1;
        has_pu_weight = false;
        has_lepton_weight = false;
    }

    /*---Run2---*/
    void CalculatePuWeight(const ntuple::Sync& eventInfo)
    {
        if(has_pu_weight)
            throw exception("PU weight is already calculated.");
        if(!apply_pu_weight)
            throw exception("PU weight should not be applied.");

        const double nPU = eventInfo.npu;
        const Int_t bin = pu_weights->FindBin(nPU);
        const Int_t maxBin = pu_weights->FindBin(max_available_pu);
        const bool goodBin = bin >= 1 && bin <= maxBin;
        PUweight = goodBin ? pu_weights->GetBinContent(bin) : default_pu_weight;

        eventWeight *= PUweight;
        has_pu_weight = true;
    }

    void CalculateLeptonWeights (const ntuple::Sync& eventInfo){

        if(has_lepton_weight)
          throw exception("Lepton weights are already calculated.");

        isoIdLeptonWeight   = isoIdScaleFactor->get_ScaleFactor(eventInfo.pt_1,eventInfo.eta_1);
        triggerLeptonWeight = triggerScaleFactor->get_ScaleFactor(eventInfo.pt_1,eventInfo.eta_1);

        eventWeight *= isoIdLeptonWeight * triggerLeptonWeight;
        has_lepton_weight = true;
    }
    /*-----------*/

    bool HasPileUpWeight()  const { return !apply_pu_weight || has_pu_weight; }
    bool HasLeptonWeights() const { return has_lepton_weight; }

    bool IsData() const { return is_data; }
    bool IsEmbedded() const { return is_embedded; }

    double GetPileUpWeight() const
    {
        if(!HasPileUpWeight())
            throw exception("Pile-up weight is not calculated yet.");
        return PUweight;
    }

    double GetIsoIdLeptonWeight() const
    {
      if(!HasLeptonWeights()) throw exception("Lepton weights are not calculated yet.");
      return isoIdLeptonWeight;
    }

    double GetTriggerLeptonWeight() const
    {
      if(!HasLeptonWeights()) throw exception("Lepton weights are not calculated yet.");
      return triggerLeptonWeight;
    }

    double GetLeptonWeights() const
    {
      if(!HasLeptonWeights()) throw exception("Lepton weights are not calculated yet.");
      return triggerLeptonWeight * isoIdLeptonWeight ;
    }

    double GetEventWeight() const { return eventWeight; }

private:
    static std::shared_ptr<TH1F> LoadPUWeights(const std::string& reweightFileName)
    {
        auto reweightFile = root_ext::OpenRootFile(reweightFileName);
        TH1F* originalWeights = root_ext::ReadObject<TH1F>(*reweightFile, "lumiWeights");
        std::shared_ptr<TH1F> weights_clone(root_ext::CloneObject(*originalWeights, "PUweights", true));
        return weights_clone;
    }

private:
    bool is_data, is_embedded, apply_pu_weight;
    double max_available_pu, default_pu_weight;
    std::shared_ptr<TH1F> pu_weights;
    bool has_pu_weight, has_lepton_weight;
    double eventWeight, PUweight, isoIdLeptonWeight, triggerLeptonWeight;
    htt_utilities::ScaleFactor * isoIdScaleFactor, *triggerScaleFactor;
};

} // namespace analysis
