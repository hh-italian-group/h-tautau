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

    JetInfo(const JetCandidate& jet);
};

struct BTagReaderInfo {
    using Reader = btag_calibration::BTagCalibrationReader;
    using ReaderPtr = std::shared_ptr<Reader>;
    using HistPtr = std::shared_ptr<TH2D>;
    using JetFlavor = btag_calibration::BTagEntry::JetFlavor;
    using FilePtr = std::shared_ptr<TFile>;

    ReaderPtr reader;
    JetFlavor flavor;
    HistPtr eff_hist;

    BTagReaderInfo(ReaderPtr _reader, JetFlavor _flavor, FilePtr file, DiscriminatorWP wp);
    void Eval(JetInfo& jetInfo, const std::string& unc_name, double btag);

private:
    double GetEfficiency(double pt, double eta) const;
};

} // namespace detail

class BTagWeight : public IWeightProvider {
public:
    using BTagCalibration = btag_calibration::BTagCalibration;
    using BTagEntry = btag_calibration::BTagEntry;
    using OperatingPoint = BTagEntry::OperatingPoint;
    using JetFlavor = BTagEntry::JetFlavor;
    using JetInfo = detail::JetInfo;
    using JetInfoVector = std::vector<JetInfo>;
    using ReaderInfo = detail::BTagReaderInfo;
    using ReaderInfoPtr = std::shared_ptr<ReaderInfo>;
    using ReaderInfoMap = std::map<DiscriminatorWP, std::map<int, ReaderInfoPtr>>;
    using Reader = ReaderInfo::Reader;
    using ReaderPtr = ReaderInfo::ReaderPtr;

    BTagWeight(const std::string& bTagEffFileName, const std::string& bjetSFFileName, const BTagger& _bTagger,
               DiscriminatorWP _default_wp);

    virtual double Get(EventInfo& event) const override;
    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override;

    double Get(EventInfo& eventInfo, DiscriminatorWP wp, UncertaintySource unc_source,
               UncertaintyScale unc_scale) const;
    std::map<UncertaintyScale, std::vector<float>> GetEvtWeightShifted(EventInfo& eventInfo, DiscriminatorWP wp,
                                                                       UncertaintySource unc_source,
                                                                       UncertaintyScale unc_scale,
                                                                       bool apply_JES = true) const;
private:
    static std::string GetUncertantyName(UncertaintyScale unc);
    static double GetBtagWeight(const JetInfoVector& jetInfos);
    ReaderInfo& GetReader(DiscriminatorWP wp, int hadronFlavour) const;

private:
    BTagCalibration calib;
    ReaderInfoMap readerInfos;
    BTagger bTagger;
    DiscriminatorWP default_wp;
    std::vector<UncertaintySource> btag_sources;
    std::vector<std::string> sist_names; 
};

} // namespace mc_corrections
} // namespace analysis
