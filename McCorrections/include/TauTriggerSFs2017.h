#pragma once

#include <string>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

namespace analysis {
namespace mc_corrections {

class TauTriggerSFs2017
{
private:
    const TH1* loadTH1(const TFile* inputFile, const std::string& histogramName);
    const TH2* loadTH2(const TFile* inputFile, const std::string& histogramName);

public:
    enum { kCentral, kStatUp, kStatDown };

    TauTriggerSFs2017(const std::string& inputFileName, const std::string& inputFileNameOld,
                      const std::string& tauMVAWP = "medium");
    ~TauTriggerSFs2017();

    double getEfficiency(double pt, double eta, double phi, const TH1* effHist, const TH2* etaPhi,
                         const TH2* etaPhiAvg, int central_or_shift = TauTriggerSFs2017::kCentral);
    double getDiTauEfficiencyData(double pt, double eta, double phi, int central_or_shift);
    double getDiTauEfficiencyMC(double pt, double eta, double phi, int central_or_shift);
    double getDiTauScaleFactor(double pt, double eta, double phi, int central_or_shift);
    double getMuTauEfficiencyData(double pt, double eta, double phi, int central_or_shift);
    double getMuTauEfficiencyMC(double pt, double eta, double phi, int central_or_shift);
    double getMuTauScaleFactor(double pt, double eta, double phi, int central_or_shift);
    double getETauEfficiencyData(double pt, double eta, double phi, int central_or_shift);
    double getETauEfficiencyMC(double pt, double eta, double phi, int central_or_shift);
    double getETauScaleFactor(double pt, double eta, double phi, int central_or_shift);

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

} // namespace mc_corrections
} // namespace analysis
