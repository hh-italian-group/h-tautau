/*! b tag weight.
Original code: https://github.com/ajgilbert/ICHiggsTauTau/blob/master/Analysis/Utilities/src/BTagWeight.cc
This file is part of https://github.com/hh-italian-group/hh-bbtautau. */

#include "h-tautau/McCorrections/include/BTagWeight.h"

namespace analysis {
namespace mc_corrections {

namespace detail {

JetInfo::JetInfo(const JetCandidate& jet) :
        eff(0.), SF(0.)
{
    pt  = jet.GetMomentum().pt();
    eta = jet.GetMomentum().eta();
    hadronFlavour = jet->hadronFlavour();
}

BTagReaderInfo::BTagReaderInfo(ReaderPtr _reader, JetFlavor _flavor, FilePtr file, DiscriminatorWP wp) :
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

void BTagReaderInfo::Eval(JetInfo& jetInfo, const std::string& unc_name, double btag)
{
    jetInfo.SF  = reader->eval_auto_bounds(unc_name, flavor, static_cast<float>(jetInfo.eta),
                                           static_cast<float>(jetInfo.pt), static_cast<float>(btag));
    jetInfo.eff = GetEfficiency(jetInfo.pt, std::abs(jetInfo.eta));
}

double BTagReaderInfo::GetEfficiency(double pt, double eta) const
{
    int xBin = eff_hist->GetXaxis()->FindFixBin(pt);
    xBin = std::min(eff_hist->GetXaxis()->GetNbins(), std::max(1, xBin));
    int yBin = eff_hist->GetYaxis()->FindFixBin(eta);
    yBin = std::min(eff_hist->GetYaxis()->GetNbins(), std::max(1, yBin));
    return eff_hist->GetBinContent(xBin, yBin);
}

} // namespace detail

BTagWeight::BTagWeight(const std::string& bTagEffFileName, const std::string& bjetSFFileName, const BTagger& _bTagger,
                       DiscriminatorWP _default_wp) :
    calib(ToString(_bTagger.GetBaseTagger()), bjetSFFileName), bTagger(_bTagger), default_wp(_default_wp)
{
    static const std::map<DiscriminatorWP, OperatingPoint> op_map = {
        { DiscriminatorWP::Loose, BTagEntry::OP_LOOSE }, { DiscriminatorWP::Medium, BTagEntry::OP_MEDIUM },
        { DiscriminatorWP::Tight, BTagEntry::OP_TIGHT }, { DiscriminatorWP::VVVLoose, BTagEntry::OP_RESHAPING }
    };
    static const std::map<JetFlavor, int> jet_flavors {
        { BTagEntry::FLAV_B, 5 }, { BTagEntry::FLAV_C, 4 }, { BTagEntry::FLAV_UDSG, 0 },
    };
    static const std::vector<std::string> unc_scale_names_comb = { "up", "down" };

    //From https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration#Systematic_uncertainties
    // Owing to the iterative approach, "jes", "lf", "hf", "hfstats1/2", and "lfstats1/2" uncertainties are applied
    //to both b and udsg jets. For c-flavored jets, only "cferr1/2" uncertainties are applied.
    btag_sources = { UncertaintySource::btag_lf, UncertaintySource::btag_hf, UncertaintySource::btag_hfstats1,
                     UncertaintySource::btag_hfstats2, UncertaintySource::btag_lfstats1, UncertaintySource::btag_lfstats2,
                     UncertaintySource::btag_cferr1, UncertaintySource::btag_cferr2, UncertaintySource::JetFull_Total,
                     UncertaintySource::JetReduced_Total };

     // Systematics names for iterative btag weights
    syst_names = { "up_jes", "up_lf", "up_hf", "up_hfstats1", "up_hfstats2", "up_lfstats1", "up_lfstats2",
                   "up_cferr1", "up_cferr2", "down_jes", "down_lf", "down_hf", "down_hfstats1", "down_hfstats2",
                   "down_lfstats1", "down_lfstats2", "down_cferr1", "down_cferr2"
    };

    if(!op_map.count(default_wp))
        throw exception("BTagWeight: default working point %1% is not supported.") % default_wp;

    auto bTagEffFile = root_ext::OpenRootFile(bTagEffFileName);

    for(const auto& [wp, op] : op_map) {
        for(const auto& [flavor, flavor_id] : jet_flavors) {
            const auto& unc_scale_names = op == BTagEntry::OP_RESHAPING ? syst_names
                                                                        : unc_scale_names_comb;
            std::string measurementType = flavor_id == 0 ? "incl" : "comb";
            if(op == BTagEntry::OP_RESHAPING)
                measurementType = "iterativefit";
            auto reader = std::make_shared<Reader>(op, "central", unc_scale_names);
            reader->load(calib, flavor, measurementType);
            const DiscriminatorWP eff_wp = op == BTagEntry::OP_RESHAPING ? DiscriminatorWP::Medium : wp;
            readerInfos[wp][flavor_id] = std::make_shared<ReaderInfo>(reader, flavor, bTagEffFile, eff_wp);
        }
    }
}

double BTagWeight::Get(EventInfo& eventInfo) const
{
    return Get(eventInfo, default_wp, false, UncertaintySource::None, UncertaintyScale::Central, true);
}

double BTagWeight::Get(const ntuple::ExpressEvent& /*event*/) const
{
    throw exception("ExpressEvent is not supported in BTagWeight::Get.");
}

double BTagWeight::Get(EventInfo& eventInfo, DiscriminatorWP wp, bool use_iterative_fit, UncertaintySource unc_source,
                       UncertaintyScale unc_scale,  bool apply_JES) const
{
    static const std::map<std::pair<UncertaintyScale,UncertaintySource>, std::string> iter_unc_scales = {
        { {UncertaintyScale::Up, UncertaintySource::btag_lf}, "up_lf" },
        { {UncertaintyScale::Down, UncertaintySource::btag_lf}, "down_lf" },
        { {UncertaintyScale::Up, UncertaintySource::btag_hf}, "up_hf" },
        { {UncertaintyScale::Down, UncertaintySource::btag_hf}, "down_hf" },
        { {UncertaintyScale::Up, UncertaintySource::btag_hfstats1}, "up_hfstats1" },
        { {UncertaintyScale::Down, UncertaintySource::btag_hfstats1}, "down_hfstats1" },
        { {UncertaintyScale::Up, UncertaintySource::btag_hfstats2}, "up_hfstats2" },
        { {UncertaintyScale::Down, UncertaintySource::btag_hfstats2}, "down_hfstats2" },
        { {UncertaintyScale::Up, UncertaintySource::btag_lfstats1}, "up_lfstats1" },
        { {UncertaintyScale::Down, UncertaintySource::btag_lfstats1}, "down_lfstats1" },
        { {UncertaintyScale::Up, UncertaintySource::btag_lfstats2}, "up_lfstats2" },
        { {UncertaintyScale::Down, UncertaintySource::btag_lfstats2}, "down_lfstats2" },
        { {UncertaintyScale::Up, UncertaintySource::btag_cferr1}, "up_cferr1" },
        { {UncertaintyScale::Down, UncertaintySource::btag_cferr1}, "down_cferr1" },
        { {UncertaintyScale::Up, UncertaintySource::btag_cferr2}, "up_cferr2" },
        { {UncertaintyScale::Down, UncertaintySource::btag_cferr2}, "down_cferr2" },
        { {UncertaintyScale::Up, UncertaintySource::JetFull_Total}, "up_jes" },
        { {UncertaintyScale::Down, UncertaintySource::JetFull_Total}, "down_jes" },
        { {UncertaintyScale::Up, UncertaintySource::JetReduced_Total}, "up_jes" },
        { {UncertaintyScale::Down, UncertaintySource::JetReduced_Total}, "down_jes" },
};

   UncertaintySource source = UncertaintySource::None;
   UncertaintyScale scale = UncertaintyScale::Central;
   if(unc_source == UncertaintySource::Eff_b)
       scale = unc_scale;
   const std::string unc_name = GetUncertantyName(scale);

   const DiscriminatorWP reader_wp = use_iterative_fit ? DiscriminatorWP::VVVLoose : wp;

   JetInfoVector jetInfos;
   double SF = 1.;
   for(const auto& jet : eventInfo.GetCentralJets()) {
       JetInfo jetInfo(*jet);
       if(!(jetInfo.pt > bTagger.PtCut() && std::abs(jetInfo.eta) < bTagger.EtaCut())) continue;
       if(use_iterative_fit){
           //For c-flavored jets, only "cferr1/2" uncertainties are applied.
           if(jetInfo.hadronFlavour == 4 && (unc_source == UncertaintySource::btag_cferr1 ||
               unc_source == UncertaintySource::btag_cferr2)) {
               scale = unc_scale;
               source = unc_source;
           }
           else if(std::count(btag_sources.begin(), btag_sources.end(), unc_source) ||
               (apply_JES && (unc_source == UncertaintySource::JetFull_Total ||
               unc_source == UncertaintySource::JetReduced_Total))) {
                   scale = unc_scale;
                   source = unc_source;
           }
       }
       const std::string unc_name_string = (use_iterative_fit && source != UncertaintySource::None) ?
           iter_unc_scales.at(std::make_pair(scale, source)) : unc_name;
       GetReader(reader_wp, jetInfo.hadronFlavour).Eval(jetInfo, unc_name_string, bTagger.BTag(**jet, true));
       SF *= jetInfo.SF;
       jetInfo.bTagOutcome = bTagger.Pass(**jet, wp);
       jetInfos.push_back(jetInfo);
   }
   return use_iterative_fit ? SF : GetBtagWeight(jetInfos);
}

std::string BTagWeight::GetUncertantyName(UncertaintyScale unc)
{
    std::string unc_name = ToString(unc);
    std::transform(unc_name.begin(), unc_name.end(), unc_name.begin(), ::tolower);
    return unc_name;
}

double BTagWeight::GetBtagWeight(const JetInfoVector& jetInfos)
{
    double MC = 1;
    double Data = 1;
    for (const auto& jetInfo : jetInfos) {
        MC *= jetInfo.bTagOutcome ? jetInfo.eff : 1 - jetInfo.eff;
        Data *= jetInfo.bTagOutcome ? jetInfo.eff * jetInfo.SF : 1 - jetInfo.eff * jetInfo.SF;
    }
    return MC != 0 ? Data/MC : 0;
}

BTagWeight::ReaderInfo& BTagWeight::GetReader(DiscriminatorWP wp, int hadronFlavour) const
{
    static const int default_flavour = 0;

    auto wp_iter = readerInfos.find(wp);
    if(wp_iter == readerInfos.end())
        throw exception("BTagWeight: working point %1% is not supported.") % wp;
    int flavour = std::abs(hadronFlavour);
    auto flav_iter = wp_iter->second.find(flavour);
    return flav_iter == wp_iter->second.end() ? *wp_iter->second.at(default_flavour) : *flav_iter->second;
}

} // namespace mc_corrections
} // namespace analysis
