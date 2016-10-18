/*! b tag weight.
Original code: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/BTagWeight.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <TH2.h>

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "BTagCalibrationStandalone.h"

namespace analysis {
namespace mc_corrections {

namespace detail {

struct JetInfo {
    double pt, eta, CSV;
    int   partonFlavour;
    double eff, SF;
    bool bTagOutcome;

    JetInfo(const ntuple::Event& event, size_t jetIndex) :
        eff(0.), SF(0.)
    {
        pt  = event.jets_p4.at(jetIndex).pt();
        eta = event.jets_p4.at(jetIndex).eta();
        CSV = event.jets_csv.at(jetIndex);
        partonFlavour = event.jets_partonFlavour.at(jetIndex);
    }
};

struct BTagReaderInfo {
    using ReaderPtr = std::shared_ptr<btag_calibration::BTagCalibrationReader>;
    using HistPtr = std::shared_ptr<TH2F>;
    using JetFlavor = btag_calibration::BTagEntry::JetFlavor;
    using FilePtr = std::shared_ptr<TFile>;

    ReaderPtr reader;
    JetFlavor flavor;
    HistPtr eff_hist;

    BTagReaderInfo(ReaderPtr _reader, JetFlavor _flavor, FilePtr file) :
        reader(_reader), flavor(_flavor)
    {
        static const std::map<JetFlavor, std::string> hist_names = {
            { btag_calibration::BTagEntry::FLAV_B, "eff_b" },
            { btag_calibration::BTagEntry::FLAV_C, "eff_c" },
            { btag_calibration::BTagEntry::FLAV_UDSG, "eff_l" },
        };

        if(!hist_names.count(flavor))
            throw exception("Jet flavor %1% not supported.") % flavor;
        eff_hist = HistPtr(root_ext::ReadCloneObject<TH2F>(*file, hist_names.at(flavor), "", true));
    }

    void Eval(JetInfo& jetInfo)
    {
        jetInfo.SF  = reader->eval(flavor, jetInfo.eta, jetInfo.pt);
        jetInfo.eff = GetEfficiency(jetInfo.pt, std::abs(jetInfo.eta));
    }

private:
    double GetEfficiency(double pt, double eta) const
    {
        const int xBin = eff_hist->GetXaxis()->FindFixBin(pt);
        const int yBin = eff_hist->GetYaxis()->FindFixBin(eta);
        return eff_hist->GetBinContent(xBin, yBin);
    }
};

} // namespace detail

class BTagWeight {
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

    BTagWeight(const std::string& bTagEffFileName, const std::string& bjetSFFileName, DiscriminatorWP wp) :
        calib("CSVv2", bjetSFFileName)
    {
        static const std::map<DiscriminatorWP, OperatingPoint> op_map = {
            { DiscriminatorWP::Loose, BTagEntry::OP_LOOSE }, { DiscriminatorWP::Medium, BTagEntry::OP_MEDIUM },
            { DiscriminatorWP::Tight, BTagEntry::OP_TIGHT }
        };
        static const std::map<JetFlavor, int> jet_flavors {
            { BTagEntry::FLAV_B, 5 }, { BTagEntry::FLAV_C, 4 }, { BTagEntry::FLAV_UDSG, 0 },
        };

        if(!op_map.count(wp))
            throw exception("Working point %1% is not supported.") % wp;

        auto bTagEffFile = root_ext::OpenRootFile(bTagEffFileName);

        ReaderPtr reader(new BTagCalibrationReader(&calib, op_map.at(wp), "mujets", "central"));
        readerInfos[jet_flavors.at(BTagEntry::FLAV_B)] =
                ReaderInfoPtr(new ReaderInfo(reader, BTagEntry::FLAV_B, bTagEffFile));
        readerInfos[jet_flavors.at(BTagEntry::FLAV_C)] =
                ReaderInfoPtr(new ReaderInfo(reader, BTagEntry::FLAV_C, bTagEffFile));

        ReaderPtr reader_light(new BTagCalibrationReader(&calib, op_map.at(wp), "incl", "central"));
        readerInfos[jet_flavors.at(BTagEntry::FLAV_UDSG)] =
                ReaderInfoPtr(new ReaderInfo(reader_light, BTagEntry::FLAV_UDSG, bTagEffFile));
    }

    double Compute(const ntuple::Event& event, double csv_cut)
    {
        JetInfoVector jetInfos;
        for (size_t jetIndex = 0; jetIndex < event.jets_p4.size(); ++jetIndex) {
            JetInfo jetInfo(event, jetIndex);
            GetReader(jetInfo.partonFlavour).Eval(jetInfo);
            jetInfo.bTagOutcome = jetInfo.CSV > csv_cut;
            jetInfos.push_back(jetInfo);
        }

        return GetBtagWeight(jetInfos);
    }

private:
    static double GetBtagWeight(const JetInfoVector& jetInfos)
    {
        double MC = 1;
        double Data = 1;
        for (const auto& jetInfo : jetInfos) {
            MC *= jetInfo.bTagOutcome ? jetInfo.eff : 1 - jetInfo.eff;
            Data *= jetInfo.bTagOutcome ? jetInfo.eff * jetInfo.SF : 1 - jetInfo.eff * jetInfo.SF;
        }
        return MC ? Data/MC : 0;
    }

    ReaderInfo& GetReader(int partonFlavour) const
    {
        static const int default_flavour = 0;
        int flavour = std::abs(partonFlavour);
        if(!readerInfos.count(flavour))
            flavour = default_flavour;
        return *readerInfos.at(flavour);
    }

private:
    BTagCalibration calib;
    ReaderInfoMap readerInfos;
};

} // namespace mc_corrections
} // namespace analysis
