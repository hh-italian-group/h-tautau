/*! b tag weight.
Original code: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/BTagWeight.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include <TH2.h>

#include "AnalysisTools/Core/include/RootExt.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "BTagCalibrationStandalone.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "WeightProvider.h"
#include "TextIO.h"

namespace analysis {
namespace mc_corrections {

namespace detail {

struct JetInfo {
    double pt, eta, CSV;
    int hadronFlavour;
    double eff, SF;
    bool bTagOutcome;

    JetInfo(const ntuple::Event& event, size_t jetIndex) :
        eff(0.), SF(0.)
    {
        pt  = event.jets_p4.at(jetIndex).pt();
        eta = event.jets_p4.at(jetIndex).eta();
        CSV = event.jets_csv.at(jetIndex);
        hadronFlavour = event.jets_hadronFlavour.at(jetIndex);
    }
};

struct BTagReaderInfo {
    using ReaderPtr = std::shared_ptr<btag_calibration::BTagCalibrationReader>;
    using HistPtr = std::shared_ptr<TH2D>;
    using JetFlavor = btag_calibration::BTagEntry::JetFlavor;
    using FilePtr = std::shared_ptr<TFile>;

    ReaderPtr reader;
    JetFlavor flavor;
    HistPtr eff_hist;

    BTagReaderInfo(ReaderPtr _reader, JetFlavor _flavor, FilePtr file, DiscriminatorWP wp) :
        reader(_reader), flavor(_flavor)
    {
        static const std::map<JetFlavor, std::string> flavor_prefixes = {
            { btag_calibration::BTagEntry::FLAV_B, "eff_b" },
            { btag_calibration::BTagEntry::FLAV_C, "eff_c" },
            { btag_calibration::BTagEntry::FLAV_UDSG, "eff_udsg" },
        };

        static const std::map<DiscriminatorWP, std::string> wp_prefixes = {
            { DiscriminatorWP::Loose, "L" }, { DiscriminatorWP::Medium, "M" },
            { DiscriminatorWP::Tight, "T" }
        };

        if(!flavor_prefixes.count(flavor))
            throw exception("Jet flavor %1% not supported.") % flavor;
        if(!wp_prefixes.count(wp))
            throw exception("B tag working point %1% not supported.") % wp;

        const std::string name = boost::str(boost::format("All/Efficiency/%1%_%2%_all")
                                            % flavor_prefixes.at(flavor) % wp_prefixes.at(wp));
        eff_hist = HistPtr(root_ext::ReadCloneObject<TH2D>(*file, name, "", true));
    }

    void Eval(JetInfo& jetInfo, const std::string& unc_name)
    {
        jetInfo.SF  = reader->eval_auto_bounds(unc_name, flavor, static_cast<float>(jetInfo.eta),
                                               static_cast<float>(jetInfo.pt));
        jetInfo.eff = GetEfficiency(jetInfo.pt, std::abs(jetInfo.eta));
    }

private:
    double GetEfficiency(double pt, double eta) const
    {
        int xBin = eff_hist->GetXaxis()->FindFixBin(pt);
        xBin = std::min(eff_hist->GetXaxis()->GetNbins(), std::max(1, xBin));
        int yBin = eff_hist->GetYaxis()->FindFixBin(eta);
        yBin = std::min(eff_hist->GetYaxis()->GetNbins(), std::max(1, yBin));
        return eff_hist->GetBinContent(xBin, yBin);
    }
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

    BTagWeight(const std::string& bTagEffFileName, const std::string& bjetSFFileName, DiscriminatorWP wp) :
        calib("CSVv2", bjetSFFileName)
    {
        static const std::map<DiscriminatorWP, std::pair<OperatingPoint, double>> op_map = {
            { DiscriminatorWP::Loose, { BTagEntry::OP_LOOSE, cuts::btag_2016::CSVv2L } },
            { DiscriminatorWP::Medium, { BTagEntry::OP_MEDIUM, cuts::btag_2016::CSVv2M } },
            { DiscriminatorWP::Tight, { BTagEntry::OP_TIGHT, cuts::btag_2016::CSVv2T } }
        };
        static const std::map<JetFlavor, int> jet_flavors {
            { BTagEntry::FLAV_B, 5 }, { BTagEntry::FLAV_C, 4 }, { BTagEntry::FLAV_UDSG, 0 },
        };

        if(!op_map.count(wp))
            throw exception("Working point %1% is not supported.") % wp;
        csv_cut = op_map.at(wp).second;
        auto bTagEffFile = root_ext::OpenRootFile(bTagEffFileName);

        ReaderPtr reader_b(new BTagCalibrationReader(op_map.at(wp).first, "central", {"up", "down"}));
        reader_b->load(calib, BTagEntry::FLAV_B, "comb");
        readerInfos[jet_flavors.at(BTagEntry::FLAV_B)] =
                ReaderInfoPtr(new ReaderInfo(reader_b, BTagEntry::FLAV_B, bTagEffFile, wp));

        ReaderPtr reader_c(new BTagCalibrationReader(op_map.at(wp).first, "central", {"up", "down"}));
        reader_c->load(calib, BTagEntry::FLAV_C, "comb");
        readerInfos[jet_flavors.at(BTagEntry::FLAV_C)] =
                ReaderInfoPtr(new ReaderInfo(reader_c, BTagEntry::FLAV_C, bTagEffFile, wp));

        ReaderPtr reader_light(new BTagCalibrationReader(op_map.at(wp).first, "central", {"up", "down"}));
        reader_light->load(calib, BTagEntry::FLAV_UDSG, "incl");
        readerInfos[jet_flavors.at(BTagEntry::FLAV_UDSG)] =
                ReaderInfoPtr(new ReaderInfo(reader_light, BTagEntry::FLAV_UDSG, bTagEffFile, wp));
    }

    virtual double Get(const ntuple::Event& event) const override
    {
        return GetEx(event, UncertaintyScale::Central);
    }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in BTagWeight::Get.");
    }

    double GetEx(const ntuple::Event& event, UncertaintyScale unc) const
    {
        const std::string unc_name = GetUncertantyName(unc);

        JetInfoVector jetInfos;
        for (size_t jetIndex = 0; jetIndex < event.jets_p4.size(); ++jetIndex) {
            JetInfo jetInfo(event, jetIndex);
            if(std::abs(jetInfo.eta) >= cuts::btag_2016::eta) continue;
            GetReader(jetInfo.hadronFlavour).Eval(jetInfo, unc_name);
            jetInfo.bTagOutcome = std::abs(jetInfo.eta) < cuts::btag_2016::eta && jetInfo.CSV > csv_cut;
            jetInfos.push_back(jetInfo);
        }

        return GetBtagWeight(jetInfos);
    }

private:

    static std::string GetUncertantyName(UncertaintyScale unc)
    {
        std::string unc_name = ToString(unc);
        std::transform(unc_name.begin(), unc_name.end(), unc_name.begin(), ::tolower);
        return unc_name;
    }

    static double GetBtagWeight(const JetInfoVector& jetInfos)
    {
        double MC = 1;
        double Data = 1;
        for (const auto& jetInfo : jetInfos) {
            MC *= jetInfo.bTagOutcome ? jetInfo.eff : 1 - jetInfo.eff;
            Data *= jetInfo.bTagOutcome ? jetInfo.eff * jetInfo.SF : 1 - jetInfo.eff * jetInfo.SF;
        }
        return MC != 0 ? Data/MC : 0;
    }

    ReaderInfo& GetReader(int hadronFlavour) const
    {
        static const int default_flavour = 0;
        int flavour = std::abs(hadronFlavour);
        if(!readerInfos.count(flavour))
            flavour = default_flavour;
        return *readerInfos.at(flavour);
    }

private:
    BTagCalibration calib;
    ReaderInfoMap readerInfos;
    double csv_cut;
};

} // namespace mc_corrections
} // namespace analysis
