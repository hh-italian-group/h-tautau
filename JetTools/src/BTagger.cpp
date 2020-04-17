/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Cuts/include/btag_Run2.h"

namespace analysis {

BTagger::BTagger(Period _period, BTaggerKind _tagger) :
    period(_period), tagger(_tagger)
{
    static const std::map<Period, std::map<BTaggerKind, std::map<DiscriminatorWP, double>>> working_points = {
        {Period::Run2018, {
            {BTaggerKind::DeepCSV, {
                {DiscriminatorWP::Loose, cuts::btag_Run2::Run2018::DeepCSVL},
                {DiscriminatorWP::Medium, cuts::btag_Run2::Run2018::DeepCSVM},
                {DiscriminatorWP::Tight, cuts::btag_Run2::Run2018::DeepCSVT}
            }},
            {BTaggerKind::DeepFlavour, {
                {DiscriminatorWP::Loose, cuts::btag_Run2::Run2018::deepFlavourL},
                {DiscriminatorWP::Medium, cuts::btag_Run2::Run2018::deepFlavourM},
                {DiscriminatorWP::Tight, cuts::btag_Run2::Run2018::deepFlavourT}
            }},
            {BTaggerKind::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_Run2::pt}
            }}
        }},
        {Period::Run2017, {
            {BTaggerKind::DeepCSV, {
                {DiscriminatorWP::Loose, cuts::btag_Run2::Run2017::deepCSVv2L},
                {DiscriminatorWP::Medium, cuts::btag_Run2::Run2017::deepCSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_Run2::Run2017::deepCSVv2T}
            }},
            {BTaggerKind::CSV, {
                {DiscriminatorWP::Loose, cuts::btag_Run2::Run2017::CSVv2L},
                {DiscriminatorWP::Medium, cuts::btag_Run2::Run2017::CSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_Run2::Run2017::CSVv2T}
            }},
            {BTaggerKind::DeepFlavour, {
                {DiscriminatorWP::Loose, cuts::btag_Run2::Run2017::deepFlavourL},
                {DiscriminatorWP::Medium, cuts::btag_Run2::Run2017::deepFlavourM},
                {DiscriminatorWP::Tight, cuts::btag_Run2::Run2017::deepFlavourT}
            }},
            {BTaggerKind::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_Run2::pt}
            }}
        }},
        {Period::Run2016, {
            {BTaggerKind::CSV, {
                {DiscriminatorWP::Loose, cuts::btag_Run2::Run2016::CSVv2L},
                {DiscriminatorWP::Medium, cuts::btag_Run2::Run2016::CSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_Run2::Run2016::CSVv2T}
            }},
             {BTaggerKind::DeepCSV, {
                 {DiscriminatorWP::Loose, cuts::btag_Run2::Run2016::DeepCSVL},
                 {DiscriminatorWP::Medium, cuts::btag_Run2::Run2016::DeepCSVM},
                 {DiscriminatorWP::Tight, cuts::btag_Run2::Run2016::DeepCSVT}
             }},
             {BTaggerKind::DeepFlavour, {
                 {DiscriminatorWP::Loose, cuts::btag_Run2::Run2016::deepFlavourL},
                 {DiscriminatorWP::Medium, cuts::btag_Run2::Run2016::deepFlavourM},
                 {DiscriminatorWP::Tight, cuts::btag_Run2::Run2016::deepFlavourT}
             }},
            {BTaggerKind::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_Run2::pt}
            }}

        }}
    };
    base_tagger = tagger;
    if(tagger == BTaggerKind::HHbtag) base_tagger = BTaggerKind::DeepFlavour;

    auto period_iter = working_points.find(period);
    if(period_iter == working_points.end())
        throw exception("Period %1% is not supported.") % period;
    auto wp_iter = period_iter->second.find(base_tagger);
    if(wp_iter == period_iter->second.end())
        throw exception("Jet tagger %1% is not supported.") % base_tagger;
    cut = &wp_iter->second;
}

double BTagger::BTag(const ntuple::TupleJet& jet, bool use_base_tagger) const
{
    const auto sel_tagger = use_base_tagger ? base_tagger : tagger;
    if(sel_tagger == BTaggerKind::Pt) return jet.p4().Pt();
    else if(sel_tagger == BTaggerKind::DeepCSV) return jet.deepcsv();
    else if (sel_tagger == BTaggerKind::CSV) return jet.csv();
    else if(sel_tagger == BTaggerKind::DeepFlavour) return jet.deepFlavour();
    else if(sel_tagger == BTaggerKind::HHbtag) return jet.hh_btag();
    throw exception("BTagger::BTag: %1% tagger is not supported") % sel_tagger;
}

bool BTagger::Pass(const ntuple::TupleJet& jet, DiscriminatorWP wp) const
{
    if(!cut->count(wp))
        throw exception("Working point %1% is not supported.") % wp;
    return BTag(jet, true) > cut->at(wp);
}

double BTagger::PtCut() const { return cuts::btag_Run2::pt; }
double BTagger::EtaCut() const { return cuts::btag_Run2::eta; }

Period BTagger::GetPeriod() const { return period; }
BTaggerKind BTagger::GetTagger() const { return tagger; }
BTaggerKind BTagger::GetBaseTagger() const { return base_tagger; }

}
