#pragma once

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <string>

namespace analysis {
namespace mc_corrections {

class TauTriggerSFs2017
{
private:
    const TH1* loadTH1(const TFile* inputFile, const std::string& histogramName)
    {
      const TH1* histogram = dynamic_cast<TH1*>((const_cast<TFile*>(inputFile))->Get(histogramName.data()));
      if ( !histogram ) {
        std::cerr << "Failed to load histogram = '" << histogramName << "' from input file !!" << std::endl;
        assert(0);
      }
      return histogram;
    }

    const TH2* loadTH2(const TFile* inputFile, const std::string& histogramName)
    {
      const TH2* histogram = dynamic_cast<TH2*>((const_cast<TFile*>(inputFile))->Get(histogramName.data()));
      if ( !histogram ) {
        std::cerr << "Failed to load histogram = '" << histogramName << "' from input file !!" << std::endl;
        assert(0);
      }
      return histogram;
    }

public:


  enum { kCentral, kStatUp, kStatDown };


  TauTriggerSFs2017(const std::string& inputFileName, const std::string& inputFileNameOld, const std::string& tauMVAWP = "medium")
  : inputFileName_(inputFileName), inputFileNameOld_(inputFileNameOld), tauMVAWP_(tauMVAWP)
    {
        inputFile_ = new TFile(inputFileName_.data());
        if ( !inputFile_ ) {
            std::cerr << "Failed to open input file = '" << inputFileName_ << "' !!" << std::endl;
            assert(0);
        }

        inputFileOld_ = new TFile(inputFileNameOld_.data());
          if ( !inputFileOld_ ) {
            std::cerr << "Failed to open input file = '" << inputFileNameOld_ << "' !!" << std::endl;
            assert(0);
          }

        // load the TH1s containing the bin by bin values
        diTauData_ = loadTH1(inputFile_, Form("hist_diTauTriggerEfficiency_%sTauMVA_DATA", tauMVAWP_.data()));
        diTauMC_ = loadTH1(inputFile_, Form("hist_diTauTriggerEfficiency_%sTauMVA_MC", tauMVAWP_.data()));
        eTauData_ = loadTH1(inputFile_, Form("hist_ETauTriggerEfficiency_%sTauMVA_DATA", tauMVAWP_.data()));
        eTauMC_ = loadTH1(inputFile_, Form("hist_ETauTriggerEfficiency_%sTauMVA_MC", tauMVAWP_.data()));
        muTauData_ = loadTH1(inputFile_, Form("hist_MuTauTriggerEfficiency_%sTauMVA_DATA", tauMVAWP_.data()));
        muTauMC_ = loadTH1(inputFile_, Form("hist_MuTauTriggerEfficiency_%sTauMVA_MC", tauMVAWP_.data()));

        // FIXME: Use the eta-phi efficiency corrections from pre-re-miniaod branch
        // Only medium, tight, and vtight are provided and they are from MVA ID
        std::string tmpMVAWP = tauMVAWP_.data();
        if (tmpMVAWP == "vvloose" || tmpMVAWP == "vloose" || tmpMVAWP == "loose") tmpMVAWP = "medium";
        if (tmpMVAWP == "vvtight")  tmpMVAWP = "vtight";

        // load the TH2s containing the eta phi efficiency corrections
        diTauEtaPhiData_ = loadTH2(inputFileOld_, Form("diTau_%s_DATA", tmpMVAWP.data()));
        diTauEtaPhiMC_ = loadTH2(inputFileOld_, Form("diTau_%s_MC", tmpMVAWP.data()));
        eTauEtaPhiData_ = loadTH2(inputFileOld_, Form("eTau_%s_DATA", tmpMVAWP.data()));
        eTauEtaPhiMC_ = loadTH2(inputFileOld_, Form("eTau_%s_MC", tmpMVAWP.data()));
        muTauEtaPhiData_ = loadTH2(inputFileOld_, Form("muTau_%s_DATA", tmpMVAWP.data()));
        muTauEtaPhiMC_ = loadTH2(inputFileOld_, Form("muTau_%s_MC", tmpMVAWP.data()));

        // Eta Phi Avg
        diTauEtaPhiAvgData_ = loadTH2(inputFileOld_, Form("diTau_%s_AVG_DATA", tmpMVAWP.data()));
        diTauEtaPhiAvgMC_ = loadTH2(inputFileOld_, Form("diTau_%s_AVG_MC", tmpMVAWP.data()));
        eTauEtaPhiAvgData_ = loadTH2(inputFileOld_, Form("eTau_%s_AVG_DATA", tmpMVAWP.data()));
        eTauEtaPhiAvgMC_ = loadTH2(inputFileOld_, Form("eTau_%s_AVG_MC", tmpMVAWP.data()));
        muTauEtaPhiAvgData_ = loadTH2(inputFileOld_, Form("muTau_%s_AVG_DATA", tmpMVAWP.data()));
        muTauEtaPhiAvgMC_ = loadTH2(inputFileOld_, Form("muTau_%s_AVG_MC", tmpMVAWP.data()));
    }
    ~TauTriggerSFs2017()
    {
      delete inputFile_;
      delete inputFileOld_;
    }

    double getEfficiency(double pt, double eta, double phi, const TH1* effHist, const TH2* etaPhi, const TH2* etaPhiAvg, int central_or_shift = TauTriggerSFs2017::kCentral)
    {
        const TAxis* effHist_xAxis = effHist->GetXaxis();
        double ptMin = effHist_xAxis->GetXmin() + 1.e-1;
        double ptMax = effHist_xAxis->GetXmax() - 1.e-1;
        double pt_checked = pt;
        if ( pt_checked > ptMax ) pt_checked = ptMax;
        if ( pt_checked < ptMin ) pt_checked = ptMin;
        int effHist_idxBin = (const_cast<TH1*>(effHist))->FindBin(pt_checked);
        assert(effHist_idxBin >= 1 && effHist_idxBin <= effHist->GetNbinsX());
        double eff = effHist->GetBinContent(effHist_idxBin);

        double effErr = effHist->GetBinError(effHist_idxBin);
        switch(central_or_shift)
        {
          case TauTriggerSFs2017::kCentral:                 break;
          case TauTriggerSFs2017::kStatUp:   eff += effErr; break;
          case TauTriggerSFs2017::kStatDown: eff -= effErr; break;
          default: assert(0);
        }

        // Adjust SF based on (eta, phi) location
        // keep eta barrel boundaries within SF region
        // but, for taus outside eta limits or with unralistic
        // phi values, return zero SF
        const TAxis* etaPhiAvg_xAxis = etaPhiAvg->GetXaxis();
        double etaMin = etaPhiAvg_xAxis->GetXmin() + 1.e-2;
        double etaMax = etaPhiAvg_xAxis->GetXmax() - 1.e-2;
        double eta_checked = eta;
        if ( eta_checked > etaMax ) eta_checked = etaMax;
        if ( eta_checked < etaMin ) eta_checked = etaMin;
        int etaPhiAvg_idxBinX = etaPhiAvg_xAxis->FindBin(eta_checked);
        assert(etaPhiAvg_idxBinX >= 1 && etaPhiAvg_idxBinX <= etaPhiAvg_xAxis->GetNbins());
        const TAxis* etaPhiAvg_yAxis = etaPhiAvg->GetYaxis();
        int etaPhiAvg_idxBinY = etaPhiAvg_yAxis->FindBin(phi);
        assert(etaPhiAvg_idxBinY >= 1 && etaPhiAvg_idxBinY <= etaPhiAvg_yAxis->GetNbins());
        double effCorr_etaPhi = etaPhi->GetBinContent((const_cast<TH2*>(etaPhi))->FindBin(eta_checked, phi));
        double effCorr_etaPhiAvg = etaPhiAvg->GetBinContent((const_cast<TH2*>(etaPhiAvg))->FindBin(eta_checked, phi));
        if ( effCorr_etaPhiAvg <= 0. ) {
          std::cerr << Form("One of the provided tau (eta, phi) values (%3.3f, %3.3f) is outside the boundary of triggering taus", eta, phi) << std::endl;
          std::cerr << "Returning efficiency = 0.0" << std::endl;
          return 0.;
        }
        eff *= (effCorr_etaPhi/effCorr_etaPhiAvg);
        if ( eff > 1. ) eff = 1.;
        return eff;
      }

    double getDiTauEfficiencyData(double pt, double eta, double phi, int central_or_shift)
    {
        return getEfficiency(pt, eta, phi, diTauData_, diTauEtaPhiData_, diTauEtaPhiAvgData_, central_or_shift);
    }

    double getDiTauEfficiencyMC(double pt, double eta, double phi, int central_or_shift)
    {
        return getEfficiency(pt, eta, phi, diTauMC_, diTauEtaPhiMC_, diTauEtaPhiAvgMC_, central_or_shift);
    }

    double getDiTauScaleFactor(double pt, double eta, double phi, int central_or_shift)
    {
        double effData = getDiTauEfficiencyData(pt, eta, phi, central_or_shift);
        double effMC = getDiTauEfficiencyMC(pt, eta, phi, central_or_shift);
        if ( effMC < 1e-5 ) {
            std::cerr << "Eff MC is suspiciously low. Please contact Tau POG." << std::endl;
            std::cerr << Form(" - DiTau Trigger SF for Tau MVA: %s   pT: %3.3f   eta: %3.3f   phi: %3.3f", tauMVAWP_.data(), pt, eta, phi) << std::endl;
            std::cerr << Form(" - MC Efficiency = %3.3f", effMC) << std::endl;
            return 0.;
        }
        double sf = (effData/effMC);
        return sf;
    }

    double getMuTauEfficiencyData(double pt, double eta, double phi, int central_or_shift)
    {
        return getEfficiency(pt, eta, phi, muTauData_, muTauEtaPhiData_, muTauEtaPhiAvgData_, central_or_shift);
    }

    double getMuTauEfficiencyMC(double pt, double eta, double phi, int central_or_shift)
    {
        return getEfficiency(pt, eta, phi, muTauMC_, muTauEtaPhiMC_, muTauEtaPhiAvgMC_, central_or_shift);
    }

    double getMuTauScaleFactor(double pt, double eta, double phi, int central_or_shift)
    {
        double effData = getMuTauEfficiencyData(pt, eta, phi, central_or_shift);
        double effMC = getMuTauEfficiencyMC(pt, eta, phi, central_or_shift);
        if ( effMC < 1e-5 ) {
            std::cerr << "Eff MC is suspiciously low. Please contact Tau POG." << std::endl;
            std::cerr << Form(" - MuTau Trigger SF for Tau MVA: %s   pT: %3.3f   eta: %3.3f   phi: %3.3f", tauMVAWP_.data(), pt, eta, phi) << std::endl;
            std::cerr << Form(" - MC Efficiency = %3.3f", effMC) << std::endl;
            return 0.;
        }
        double sf = (effData/effMC);
        return sf;
    }

    double getETauEfficiencyData(double pt, double eta, double phi, int central_or_shift)
    {
        return getEfficiency(pt, eta, phi, eTauData_, eTauEtaPhiData_, eTauEtaPhiAvgData_, central_or_shift);
    }

    double getETauEfficiencyMC(double pt, double eta, double phi, int central_or_shift)
    {
        return getEfficiency(pt, eta, phi, eTauMC_, eTauEtaPhiMC_, eTauEtaPhiAvgMC_, central_or_shift);
    }

    double getETauScaleFactor(double pt, double eta, double phi, int central_or_shift)
    {
        double effData = getETauEfficiencyData(pt, eta, phi, central_or_shift);
        double effMC = getETauEfficiencyMC(pt, eta, phi, central_or_shift);
        if ( effMC < 1e-5 ) {
            std::cerr << "Eff MC is suspiciously low. Please contact Tau POG." << std::endl;
            std::cerr << Form(" - ETau Trigger SF for Tau MVA: %s   pT: %3.3f   eta: %3.3f   phi: %3.3f", tauMVAWP_.data(), pt, eta, phi) << std::endl;
            std::cerr << Form(" - MC Efficiency = %3.3f", effMC) << std::endl;
            return 0.;
        }
        double sf = (effData/effMC);
        return sf;
    }

protected:
  std::string inputFileName_,inputFileNameOld_;
  TFile* inputFile_;
  TFile* inputFileOld_;


  std::string tauMVAWP_;

  const TH1* diTauData_;
  const TH1* diTauMC_;
  const TH1* eTauData_;
  const TH1* eTauMC_;
  const TH1* muTauData_;
  const TH1* muTauMC_;

  const TH2* diTauEtaPhiData_;
  const TH2* diTauEtaPhiMC_;
  const TH2* eTauEtaPhiData_;
  const TH2* eTauEtaPhiMC_;
  const TH2* muTauEtaPhiData_;
  const TH2* muTauEtaPhiMC_;

  const TH2* diTauEtaPhiAvgData_;
  const TH2* diTauEtaPhiAvgMC_;
  const TH2* eTauEtaPhiAvgData_;
  const TH2* eTauEtaPhiAvgMC_;
  const TH2* muTauEtaPhiAvgData_;
  const TH2* muTauEtaPhiAvgMC_;
};

}
}
