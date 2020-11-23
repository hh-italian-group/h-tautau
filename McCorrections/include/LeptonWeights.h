/* Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#pragma once

#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "TauAnalysisTools/TauTriggerSFs/interface/SFProvider.h"
#include "TauPOG/TauIDSFs/interface/TauIDSFTool.h"
#include "VBFTrigger/VBFTriggerSFs/interface/VBFTriggerSFs.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {
namespace detail {

class LeptonScaleFactors {
public:
    using SF = htt_utilities::ScaleFactor;

    LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInputSingle,
                       const std::string& triggerInputCross = "");

    double GetIdIsoSF(const LorentzVector& p4, UncertaintyScale unc_scale) const;
    double GetTriggerEff(const LorentzVector& p4, bool isData, UncertaintyScale unc_scale) const;
    double GetTriggerEffCross(const LorentzVector& p4, bool isData, UncertaintyScale unc_scale) const;
    bool HasCrossTriggers() const;

private:
    static double GetEff(const LorentzVector& p4, bool isData, UncertaintyScale unc_scale, SF& sf_provider);

private:
    std::shared_ptr<SF> idIso, triggerSingle, triggerCross;
};

class ElectronScaleFactorPOG {
public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    ElectronScaleFactorPOG(const std::string& idInput, const std::string& isoInput);

    double GetTriggerSF() const;

    template<typename LorentzVector>
    double GetIsoSF(const LorentzVector& p4) const
    {
        const Int_t bin_eta = iso_hist->GetXaxis()->FindBin(std::abs(p4.eta()));
        const Int_t bin_et = iso_hist->GetYaxis()->FindBin(p4.Et());
        double sf = iso_hist->GetBinContent(bin_eta,bin_et);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdSF(const LorentzVector& p4) const
    {
        const Int_t bin_eta = id_hist->GetXaxis()->FindBin(std::abs(p4.eta()));
        const Int_t bin_et = id_hist->GetYaxis()->FindBin(p4.Et());
        double sf = id_hist->GetBinContent(bin_eta,bin_et);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const { return GetIdSF(p4) * GetIsoSF(p4); }

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name);

private:
    HistPtr id_hist, iso_hist;
 };

 class MuonScaleFactorPOG {
 public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    MuonScaleFactorPOG(const std::string& idInput, const std::string& isoInput, const std::string& triggerInput);

    template<typename LorentzVector>
    double GetTriggerSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt = trigger_hist->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta = trigger_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf = trigger_hist->GetBinContent(bin_pt,bin_eta);

        return sf;
    }

    template<typename LorentzVector>
    double GetIsoSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt = iso_hist->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta = iso_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf = iso_hist->GetBinContent(bin_pt,bin_eta);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdSF(const LorentzVector& p4) const
    {
        const Int_t bin_pt = id_hist->GetXaxis()->FindBin(p4.pt());
        const Int_t bin_eta = id_hist->GetYaxis()->FindBin(std::abs(p4.eta()));
        double sf = id_hist->GetBinContent(bin_pt,bin_eta);

        return sf;
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const { return GetIdSF(p4) * GetIsoSF(p4); }

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name);

private:
    HistPtr id_hist, iso_hist, trigger_hist;
};

} // namespace detail

class LeptonWeights : public IWeightProvider {
public:
    using Event = ntuple::Event;

    LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_SingletriggerInput,
                  const std::string& electron_CrossTriggerInput, const std::string& muon_idIsoInput,
                  const std::string& muon_SingletriggerInput, const std::string& muon_CrossTriggerInput,
                  const std::string& tauTriggerInput, const std::string& tauVBFTriggerInput, Period period,
                  bool _is_dm_binned);

    TauIDSFTool& GetTauIdProvider(TauIdDiscriminator discr, DiscriminatorWP wp);
    const tau_trigger::SFProvider& GetTauTriggerSFProvider(Channel channel, DiscriminatorWP wp, bool is_vbf = false);
    double GetLegIdIsoWeight(LepCandidate leg, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                             DiscriminatorWP VSjet_wp, UncertaintySource unc_source, UncertaintyScale unc_scale);

    double GetIdIsoWeight(EventInfo& eventInfo, DiscriminatorWP VSe_wp, DiscriminatorWP VSmu_wp,
                          DiscriminatorWP VSjet_wp, UncertaintySource unc_source, UncertaintyScale unc_scale);

    bool ApplyIdUncertaintyScale(int decay_mode, double pt, double eta, GenLeptonMatch gen_lepton_match,
                                 LegType leg_type, UncertaintySource unc_source);

    double GetTriggerWeight(EventInfo& eventInfo, DiscriminatorWP VSjet_wp, UncertaintySource unc_source,
                            UncertaintyScale unc_scale);
    double GetTriggerPrescaleWeight(EventInfo& eventInfo) const;

    double GetTriggerEfficiency(EventInfo& eventInfo, bool isData, DiscriminatorWP VSjet_wp,
                                UncertaintySource unc_source, UncertaintyScale unc_scale, bool& same_as_central);
    double GetVBFTriggerEfficiency(EventInfo& eventInfo, bool isData,UncertaintySource unc_source,
                                   UncertaintyScale unc_scale, bool& same_as_central);
    double GetCustomTauSF(const LepCandidate& leg, UncertaintySource unc_source, UncertaintyScale unc_scale,
                          Channel channel);

    virtual double Get(EventInfo& eventInfo) const override;
    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override;

private:
    detail::LeptonScaleFactors electronSF, muonSF;
    std::string tauTriggerInput, tauVBFTriggerInput;
    Period period;
    bool is_dm_binned;
    std::map<TauIdDiscriminator, std::map<DiscriminatorWP, std::shared_ptr<TauIDSFTool>>> tau_sf_providers;
    std::map<Channel, std::map<DiscriminatorWP, std::shared_ptr<tau_trigger::SFProvider>>> tau_trigger_sf_providers;
    std::shared_ptr<VBFTriggerSFs> vbf_trigger_provider;
};

} // namespace mc_corrections
} // namespace analysis
