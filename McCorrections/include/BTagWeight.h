/*! b tag weight.
Original code: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/BTagWeight.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#pragma once

#include <TH2.h>

#include "AnalysisTools/Core/include/TextIO.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Core/include/AnalysisTypes.h"
#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/JetTools/include/BTagCalibrationStandalone.h"
#include "WeightProvider.h"

namespace analysis {
namespace mc_corrections {

namespace detail {

struct JetInfo {
    double pt, eta;
    int hadronFlavour;
    double eff, SF;
    bool bTagOutcome;

    JetInfo(EventInfo& eventInfo, size_t jetIndex);
};

struct BTagReaderInfo {
    using ReaderPtr = std::shared_ptr<btag_calibration::BTagCalibrationReader>;
    using HistPtr = std::shared_ptr<TH2D>;
    using JetFlavor = btag_calibration::BTagEntry::JetFlavor;
    using FilePtr = std::shared_ptr<TFile>;

    ReaderPtr reader;
    JetFlavor flavor;
    HistPtr eff_hist;

    BTagReaderInfo(ReaderPtr _reader, JetFlavor _flavor, FilePtr file, DiscriminatorWP wp);
    void Eval(JetInfo& jetInfo, const std::string& unc_name);

private:
    double GetEfficiency(double pt, double eta) const;
};

} // namespace detail

class BTagWeight : public IWeightProvider {
public:
    using BTagCalibration = btag_calibration::BTagCalibration;
    using BTagCalibrationReader = btag_calibration::BTagCalibrationReader;
    using BTagEntry = btag_calibration::BTagEntry;
    using OperatingPoint = BTagEntry::OperatingPoint;
    using JetFlavor = BTagEntry::JetFlavor;
    using JetInfo = detail::JetInfo;
    using JetInfoVector = std::vector<JetInfo>;
    using ReaderInfo = detail::BTagReaderInfo;
    using ReaderInfoPtr = std::shared_ptr<detail::BTagReaderInfo>;
    using ReaderInfoMap = std::map<int, ReaderInfoPtr>;
    using ReaderPtr = std::shared_ptr<btag_calibration::BTagCalibrationReader>;

    BTagWeight(const std::string& bTagEffFileName, const std::string& bjetSFFileName,Period period,
               JetOrdering ordering, DiscriminatorWP wp);

    virtual double Get(EventInfo& event) const override;
    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override;
    double GetEx(EventInfo& eventInfo, UncertaintyScale unc) const;

private:
    static std::string GetUncertantyName(UncertaintyScale unc);
    static double GetBtagWeight(const JetInfoVector& jetInfos);
    ReaderInfo& GetReader(int hadronFlavour) const;

private:
    BTagCalibration calib;
    ReaderInfoMap readerInfos;
    BTagger bTagger;
};

} // namespace mc_corrections
} // namespace analysis
