/* Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#pragma once

#include "AnalysisTools/Core/include/EventIdentifier.h"
#include "HTT-utilities/LepEffInterface/interface/ScaleFactor.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "WeightProvider.h"
#include "TauTriggerSFs2017.h"
#include "TauIdWeight.h"

namespace analysis {
namespace mc_corrections {
namespace detail {

class LeptonScaleFactors {
public:
    using SF = htt_utilities::ScaleFactor;

    LeptonScaleFactors(const std::string& idIsoInput, const std::string& triggerInputSingle, const std::string& triggerInputCross = "") :
        idIso(new SF()), triggerSingle(new SF())

    {
        idIso->init_ScaleFactor(idIsoInput);
        triggerSingle->init_ScaleFactor(triggerInputSingle);
        if(!triggerInputCross.empty()) {
            triggerCross = std::make_shared<SF>();
            triggerCross->init_ScaleFactor(triggerInputCross);
        }
    }

    template<typename LorentzVector>
    double GetIdIsoSF(const LorentzVector& p4) const
    {
        return idIso->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTriggerSF(const LorentzVector& p4 ) const
    {
        return triggerSingle->get_ScaleFactor(p4.pt(), p4.eta());
    }

    template<typename LorentzVector>
    double GetTriggerSFCross(const LorentzVector& p4 ) const
    {
        return triggerCross->get_ScaleFactor(p4.pt(), p4.eta());
    }

private:
    std::shared_ptr<SF> idIso, triggerSingle, triggerCross;
};

class ElectronScaleFactorPOG {
public:
    using Hist = TH1;
    using HistPtr = std::shared_ptr<Hist>;

    ElectronScaleFactorPOG(const std::string& idInput, const std::string& isoInput) :
        id_hist(LoadWeight(idInput,"EGamma_SF2D")),
        iso_hist(LoadWeight(isoInput,"EGamma_SF2D"))
        {}

    double GetTriggerSF() const { return 0.991; }

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

    static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name)
    {
        auto file = root_ext::OpenRootFile(weight_file_name);
        return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
    }

private:
    HistPtr id_hist, iso_hist;
 };

 class MuonScaleFactorPOG {
 public:
     using Hist = TH1;
     using HistPtr = std::shared_ptr<Hist>;

     MuonScaleFactorPOG(const std::string& idInput, const std::string& isoInput, const std::string& triggerInput) :
         id_hist(LoadWeight(idInput,"NUM_MediumID_DEN_genTracks_pt_abseta")),
         iso_hist(LoadWeight(isoInput,"NUM_TightRelIso_DEN_MediumID_pt_abseta")),
         trigger_hist(LoadWeight( triggerInput,"IsoMu27_PtEtaBins/pt_abseta_ratio"))
         {}

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

     static HistPtr LoadWeight(const std::string& weight_file_name, const std::string& hist_name)
     {
         auto file = root_ext::OpenRootFile(weight_file_name);
         return HistPtr(root_ext::ReadCloneObject<Hist>(*file, hist_name, "", true));
     }

     private:
     HistPtr id_hist, iso_hist, trigger_hist;
};

} // namespace detail

class LeptonWeights : public IWeightProvider {
public:
    using Event = ntuple::Event;

    LeptonWeights(const std::string& electron_idIsoInput, const std::string& electron_SingletriggerInput,
                  const std::string& muon_idIsoInput, const std::string& muon_SingletriggerInput ) :
        electronSF(electron_idIsoInput, electron_SingletriggerInput),
        muonSF(muon_idIsoInput, muon_SingletriggerInput)
     {
     }

    LeptonWeights(Period period, const std::string& electron_idIsoInput, const std::string& electron_SingletriggerInput,
                  const std::string& electron_CrossTriggerInput, const std::string& tauTriggerInput,
                  const std::string& muon_idIsoInput, const std::string& muon_SingletriggerInput,
                  const std::string& muon_CrossTriggerInput,
                  const std::string& electron_idInput, const std::string& electron_IsoInput,
                  const std::string& muon_idInput, const std::string& muon_IsoInput, const std::string& muon_triggerInput,
                  const std::string& tauId_input, DiscriminatorWP tau_iso_wp ) :
        electronSF(electron_idIsoInput, electron_SingletriggerInput, electron_CrossTriggerInput),
        muonSF(muon_idIsoInput, muon_SingletriggerInput, muon_CrossTriggerInput),
        tauSF(std::make_shared<TauTriggerSFs2017>(tauTriggerInput, "medium")),
        electronSFPOG(std::make_shared<detail::ElectronScaleFactorPOG>(electron_idInput, electron_IsoInput)),
        muonSFPOG(std::make_shared<detail::MuonScaleFactorPOG>(muon_idInput, muon_IsoInput, muon_triggerInput)),
        tauIdWeight(std::make_shared<TauIdWeight2017>(tauId_input, tau_iso_wp))
    {
        if(period == Period::Run2016)
            tauIdWeight = std::make_shared<TauIdWeight>(tauId_input, tau_iso_wp);
        else if(period == Period::Run2017)
            tauIdWeight = std::make_shared<TauIdWeight2017>(tauId_input, tau_iso_wp);
        else
            throw exception("Period %1% is not supported.");
    }

    double GetIdIsoWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);

        if(channel == Channel::ETau) {
            //Recommendation for Ele SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations
            double tauSF = 0;
            //Barrel ( abs(eta) < 1.460)
            if(std::abs(event.p4_1.eta()) < 1.460){
                double tauIsoMedium = 0.89;
                double AgainstEleTight = 1.80;

                tauSF = tauIsoMedium * AgainstEleTight;
            }
            // Endcaps ( abs(eta) > 1.558)
            else if(std::abs(event.p4_1.eta()) < 1.558){
                double tauIsoMedium = 0.89;
                double AgainstEleTight = 1.53;

                tauSF = tauIsoMedium * AgainstEleTight;
            }
            return (electronSF.GetIdIsoSF(event.p4_1) * tauSF);
        }

        else if(channel == Channel::MuTau) {
            //Recommendation for Muon SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017
            double tauSF = 0;
            if(0 < std::abs(event.p4_1.eta()) < 0.4){
                double tauIsoMedium = 0.89;
                double againstMuonTight = 1.17;

                tauSF = tauIsoMedium * againstMuonTight;
            }

            else if(0.4 < std::abs(event.p4_1.eta()) < 0.8){
                double tauIsoMedium = 0.89;
                double againstMuonTight = 1.29;

                tauSF = tauIsoMedium * againstMuonTight;
            }

            else if(0.8 < std::abs(event.p4_1.eta()) < 1.2){
                double tauIsoMedium = 0.89;
                double againstMuonTight = 1.14;

                tauSF = tauIsoMedium * againstMuonTight;
            }

            else if(1.2 < std::abs(event.p4_1.eta()) < 1.7){
                double tauIsoMedium = 0.89;
                double againstMuonTight = 0.93;

                tauSF = tauIsoMedium * againstMuonTight;
            }

            else if(1.7 < std::abs(event.p4_1.eta()) < 2.3){
                double tauIsoMedium = 0.89;
                double againstMuonTight = 1.61;

                tauSF = tauIsoMedium * againstMuonTight;
            }
            return muonSF.GetIdIsoSF(event.p4_1) * tauSF;
        }

        else if (channel == Channel::TauTau) {
            //Recommendation for Tau SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
            //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
            double tauIsoMedium = 0.89;
            double tauSF =  2 * tauIsoMedium;

            return tauSF;
        }

        else if(channel == Channel::MuMu)
            return muonSF.GetIdIsoSF(event.p4_1) * muonSF.GetIdIsoSF(event.p4_2);

        else
            throw exception ("channel not allowed");
    }

    double GetTriggerWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::ETau) {

            if(std::abs(event.p4_2.eta()) < 2.1){

                double SF = electronSF.GetTriggerSF(event.p4_1) *
                            (1 - tauSF->getETauScaleFactor(event.p4_2.pt(), event.p4_2.eta(), event.p4_2.phi(), TauTriggerSFs2017::kCentral)) +
                            electronSF.GetTriggerSFCross(event.p4_1) * (tauSF->getETauScaleFactor(event.p4_2.pt(), event.p4_2.eta(), event.p4_2.phi(), TauTriggerSFs2017::kCentral));
                return SF;
            }
            else
                return electronSF.GetTriggerSF(event.p4_1);
        }
        if(channel == Channel::MuTau){
            if(std::abs(event.p4_2.eta()) < 2.1){

                double SF = muonSF.GetTriggerSF(event.p4_1) *
                            (1 - tauSF->getMuTauScaleFactor(event.p4_2.pt(), event.p4_2.eta(), event.p4_2.phi(), TauTriggerSFs2017::kCentral)) +
                            muonSF.GetTriggerSFCross(event.p4_1) * tauSF->getMuTauScaleFactor(event.p4_2.pt(), event.p4_2.eta(), event.p4_2.phi(), TauTriggerSFs2017::kCentral);
                return SF;
            }
            else
                return muonSF.GetTriggerSF(event.p4_1);
            }

        if(channel == Channel::TauTau){

            double SF1 = tauSF->getDiTauScaleFactor(event.p4_1.pt(), event.p4_1.eta(), event.p4_1.phi(), TauTriggerSFs2017::kCentral);
            double SF2 = tauSF->getDiTauScaleFactor(event.p4_2.pt(), event.p4_2.eta(), event.p4_2.phi(), TauTriggerSFs2017::kCentral);

            return SF1 * SF2;
        }

        if(channel == Channel::MuMu)
            return muonSF.GetTriggerSF(event.p4_1);

        else throw exception ("channel not allowed");
    }

    virtual double Get(const Event& event) const override { return GetIdIsoWeight(event) * GetTriggerWeight(event); }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
    }

private:
    detail::LeptonScaleFactors electronSF, muonSF;
    std::shared_ptr<TauTriggerSFs2017> tauSF;
    std::shared_ptr<detail::ElectronScaleFactorPOG> electronSFPOG;
    std::shared_ptr<detail::MuonScaleFactorPOG> muonSFPOG;
    std::shared_ptr<IWeightProvider> tauIdWeight;

};

} // namespace mc_corrections
} // namespace analysis
