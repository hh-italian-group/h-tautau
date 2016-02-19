/*! Definition of the class to calculate and store different event weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "SelectionResults.h"
#include "AnalysisTools/Core/include/Tools.h"
#include "Htautau_Summer13.h"
#include "MCfinalState.h"

namespace analysis {

class EventWeights {
public:
    EventWeights(bool _is_data, bool _is_embedded, bool _apply_pu_weight, bool _apply_DM_weight,
                 const std::string& pu_reweight_file_name, double _max_available_pu, double _default_pu_weight)
        : is_data(_is_data), is_embedded(_is_embedded), apply_pu_weight(_apply_pu_weight),
          apply_DM_weight(_apply_DM_weight), max_available_pu(_max_available_pu),
          default_pu_weight(_default_pu_weight)
    {
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
        embeddedWeight = 1;
        triggerWeights.assign(SelectionResults::NumberOfLegs, 1);
        IDweights.assign(SelectionResults::NumberOfLegs, 1);
        IsoWeights.assign(SelectionResults::NumberOfLegs, 1);
        DMweights.assign(SelectionResults::NumberOfLegs, 1);
        fakeWeights.assign(SelectionResults::NumberOfLegs, 1);
        has_pu_weight = false;
        has_selection_dependent_weights = false;
        has_embedded_weight = false;
        has_gen_taus = false;
        genTaus.clear();
    }

    void CalculatePuWeight(const ntuple::Event& eventInfo)
    {
        if(has_pu_weight)
            throw exception("PU weight is already calculated.");
        if(!apply_pu_weight)
            throw exception("PU weight should not be applied.");

        const size_t bxIndex = tools::find_index(eventInfo.bunchCrossing, 0);
        if(bxIndex >= eventInfo.bunchCrossing.size())
            throw std::runtime_error("in-time BX not found");

        const double nPU = eventInfo.trueNInt.at(bxIndex);
        const Int_t bin = pu_weights->FindBin(nPU);
        const Int_t maxBin = pu_weights->FindBin(max_available_pu);
        const bool goodBin = bin >= 1 && bin <= maxBin;
        PUweight = goodBin ? pu_weights->GetBinContent(bin) : default_pu_weight;

        eventWeight *= PUweight;
        has_pu_weight = true;
    }

    void CalculateSelectionDependentWeights(const SelectionResults& selection)
    {
        if(has_selection_dependent_weights)
            throw exception("Selection dependent weights are already calculated.");
        if(is_data)
            throw exception("Selection dependent weights should not be applied");

        CalculateWeight(&EventWeights::CalculateTriggerWeight, triggerWeights, selection);
        CalculateWeight(&EventWeights::CalculateIsoWeight, IsoWeights, selection);
        CalculateWeight(&EventWeights::CalculateIdWeight, IDweights, selection);
        CalculateWeight(&EventWeights::CalculateDecayModeWeight, DMweights, selection);
        CalculateWeight(&EventWeights::CalculateFakeWeight, fakeWeights, selection);

        has_selection_dependent_weights = true;
    }

    void SetEmbeddedWeight(double _embeddedWeight)
    {
        if(has_embedded_weight)
            throw exception("Embedded weight is already set.");
        if(!is_embedded)
            throw exception("Embedded weight should not be applied.");

        embeddedWeight = _embeddedWeight;
        eventWeight *= embeddedWeight;
        has_embedded_weight = true;
    }

    void SetGenTaus(const finalState::bbTauTau& mc_final_state)
    {
        if(has_gen_taus)
            throw exception("Gen taus are already set.");

        genTaus = mc_final_state.hadronic_taus;

        has_gen_taus = true;
    }

    bool HasPileUpWeight() const { return !apply_pu_weight || has_pu_weight; }
    bool HasEmbeddedWeight() const { return !is_embedded || has_embedded_weight; }
    bool HasPartialWeight() const { return HasPileUpWeight() && HasEmbeddedWeight(); }
    bool HasSelectionDependentWeights() const { return is_data || has_selection_dependent_weights; }
    bool HasFullWeight() const { return HasPartialWeight() && HasSelectionDependentWeights(); }

    bool IsData() const { return is_data; }
    bool IsEmbedded() const { return is_embedded; }

    double GetFullWeight() const
    {
        if(!HasFullWeight())
            throw exception("Full weight is not calculated yet.");
        return eventWeight;
    }

    double GetPartialWeight() const
    {
        if(!HasPartialWeight())
            throw exception("Partial weight is not calculated yet.");
        return eventWeight;
    }

    double GetPileUpWeight() const
    {
        if(!HasPileUpWeight())
            throw exception("Pile-up weight is not calculated yet.");
        return PUweight;
    }

    double GetEmbeddedWeight() const
    {
        if(!HasEmbeddedWeight())
            throw exception("Embedded weight is not calculated yet.");
        return embeddedWeight;
    }

    double GetTriggerWeight(size_t leg_id) const { return GetWeight(triggerWeights, leg_id); }
    double GetIsoWeight(size_t leg_id) const { return GetWeight(IsoWeights, leg_id); }
    double GetIdWeight(size_t leg_id) const { return GetWeight(IDweights, leg_id); }
    double GetDecayModeWeight(size_t leg_id) const { return GetWeight(DMweights, leg_id); }
    double GetFakeWeight(size_t leg_id) const { return GetWeight(fakeWeights, leg_id); }

protected:
    virtual double CalculateTriggerWeight(CandidatePtr leg) { return 1; }
    virtual double CalculateIsoWeight(CandidatePtr leg) { return 1; }
    virtual double CalculateIdWeight(CandidatePtr leg) { return 1; }

    double CalculateDecayModeWeight(CandidatePtr leg)
    {
        using namespace cuts::Htautau_Summer13::tauCorrections;

        if(leg->GetType() != Candidate::Type::Tau || !apply_DM_weight)
            return 1;

        if(!has_gen_taus)
            throw exception("Gen taus are not set.");

        const ntuple::Tau& tau_leg = leg->GetNtupleObject<ntuple::Tau>();
        if(tau_leg.decayMode == ntuple::tau_id::kOneProng0PiZero
                && FindMatchedObjects(leg->GetMomentum(), genTaus, deltaR_matchGenParticle).size())
            return DecayModeWeight;
        return 1;
    }

    virtual double CalculateFakeWeight(CandidatePtr leg) { return 1; }

private:
    static std::shared_ptr<TH1D> LoadPUWeights(const std::string& reweightFileName)
    {
        auto reweightFile = root_ext::OpenRootFile(reweightFileName);
        TH1D* originalWeights = root_ext::ReadObject<TH1D>(*reweightFile, "weights");
        std::shared_ptr<TH1D> weights_clone(root_ext::CloneObject(*originalWeights, "PUweights", true));
        return weights_clone;
    }

    static size_t GetLegId(size_t index) { return index + 1; }
    static size_t GetIndex(size_t leg_id)
    {
        if(leg_id == 0 || leg_id > SelectionResults::NumberOfLegs)
            throw exception("Invalid leg id = ") << leg_id;
        return leg_id - 1;
    }

    template<typename CalcMethod>
    void CalculateWeight(CalcMethod method, std::vector<double>& container, const SelectionResults& selection)
    {
        for(size_t n = 0; n < SelectionResults::NumberOfLegs; ++n) {
            const size_t leg_id = GetLegId(n);
            const CandidatePtr leg = selection.GetLeg(leg_id);
            container.at(n) = (this->*method)(leg);
            eventWeight *= container.at(n);
        }
    }

    double GetWeight(const std::vector<double>& container, size_t leg_id) const
    {
        if(!HasSelectionDependentWeights())
            throw exception("Selection dependent weights are not calculated yet.");
        return container.at(GetIndex(leg_id));
    }

private:
    bool is_data, is_embedded, apply_pu_weight, apply_DM_weight;
    double max_available_pu, default_pu_weight;
    std::shared_ptr<TH1D> pu_weights;
    bool has_pu_weight, has_selection_dependent_weights, has_embedded_weight;
    double eventWeight, PUweight, embeddedWeight;
    std::vector<double> triggerWeights, IDweights, IsoWeights, DMweights, fakeWeights;
    bool has_gen_taus;
    VisibleGenObjectVector genTaus;
};

} // namespace analysis
