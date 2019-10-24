/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */


#include "h-tautau/JetTools/include/BTagger.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/Cuts/include/Btag_2016.h"

namespace analysis {

BTagger::BTagger(Period _period, JetOrdering _ordering) :
    period(_period), ordering(_ordering)
{
    static const std::map<Period, std::map<JetOrdering, std::map<DiscriminatorWP, double>>> working_points = {
        {Period::Run2018, {
            {JetOrdering::DeepCSV, {
                {DiscriminatorWP::Loose, cuts::btag_2018::deepCSVv2L},{DiscriminatorWP::Medium, cuts::btag_2018::deepCSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_2018::deepCSVv2T}
            }},
            {JetOrdering::DeepFlavour, {
                {DiscriminatorWP::Loose, cuts::btag_2018::deepFlavourL},
                {DiscriminatorWP::Medium, cuts::btag_2018::deepFlavourM},
                {DiscriminatorWP::Tight, cuts::btag_2018::deepFlavourT}
            }},
            {JetOrdering::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_2018::pt}
            }}
        }},
        {Period::Run2017, {
            {JetOrdering::DeepCSV, {
                {DiscriminatorWP::Loose, cuts::btag_2017::deepCSVv2L},{DiscriminatorWP::Medium, cuts::btag_2017::deepCSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_2017::deepCSVv2T}
            }},
            {JetOrdering::CSV, {
                {DiscriminatorWP::Loose, cuts::btag_2017::CSVv2L}, {DiscriminatorWP::Medium, cuts::btag_2017::CSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_2017::CSVv2T}
            }},
            {JetOrdering::DeepFlavour, {
                {DiscriminatorWP::Loose, cuts::btag_2017::deepFlavourL},
                {DiscriminatorWP::Medium, cuts::btag_2017::deepFlavourM},
                {DiscriminatorWP::Tight, cuts::btag_2017::deepFlavourT}
            }},
            {JetOrdering::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_2017::pt}
            }}
        }},
        {Period::Run2016, {
            {JetOrdering::CSV, {
                {DiscriminatorWP::Loose, cuts::btag_2016::CSVv2L},
                {DiscriminatorWP::Medium, cuts::btag_2016::CSVv2M},
                {DiscriminatorWP::Tight, cuts::btag_2016::CSVv2T}
            }},
             {JetOrdering::DeepCSV, {
                 {DiscriminatorWP::Loose, cuts::btag_2016::deepCSVv2L},{DiscriminatorWP::Medium, cuts::btag_2016::deepCSVv2M},
                 {DiscriminatorWP::Tight, cuts::btag_2016::deepCSVv2T}
             }},
             {JetOrdering::DeepFlavour, {
                 {DiscriminatorWP::Loose, cuts::btag_2016::deepFlavourL},
                 {DiscriminatorWP::Medium, cuts::btag_2016::deepFlavourM},
                 {DiscriminatorWP::Tight, cuts::btag_2016::deepFlavourT}
             }},
            {JetOrdering::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_2016::pt}
            }}

        }}
    };

    if(!working_points.count(period))
        throw exception("Period %1% is not supported.") % period;
    if(!working_points.at(period).count(ordering))
            throw exception("Jet Ordering %1% is not supported.") % ordering;

    cut = &working_points.at(period).at(ordering);
}

double BTagger::BTag(const ntuple::Event& event, size_t jet_index) const
{
    if(ordering==JetOrdering::Pt) return event.jets_p4.at(jet_index).Pt();
    else if(ordering==JetOrdering::DeepCSV) return event.jets_deepCsv_BvsAll.at(jet_index);
    else if (ordering==JetOrdering::CSV) return event.jets_csv.at(jet_index);
    else if(ordering==JetOrdering::DeepFlavour) return (event.jets_deepFlavour_b.at(jet_index) +
            event.jets_deepFlavour_bb.at(jet_index) + event.jets_deepFlavour_lepb.at(jet_index));
    else
        throw exception("Jet Ordering %1% is not supported") % ordering;
}

double BTagger::BTag(const ntuple::TupleJet& jet) const
{
    if(ordering==JetOrdering::Pt) return jet.p4().Pt();
    else if(ordering==JetOrdering::DeepCSV) return jet.deepcsv();
    else if(ordering==JetOrdering::CSV) return jet.csv();
    else if(ordering==JetOrdering::DeepFlavour) return jet.deepFlavour();
    else
        throw exception("Jet Ordering %1% is not supported") % ordering;
}

bool BTagger::Pass(const ntuple::Event& event, size_t jet_index, DiscriminatorWP wp) const
{
    if(!cut->count(wp))
        throw exception("Working point %1% is not supported.") %wp;
    else return BTag(event, jet_index) > cut->at(wp);
}

bool BTagger::Pass(const ntuple::TupleJet& jet, DiscriminatorWP wp) const
{
    if(!cut->count(wp))
        throw exception("Working point %1% is not supported.") %wp;
    else return BTag(jet) > cut->at(wp);
}

double BTagger::PtCut() const
{
    return period == analysis::Period::Run2017 ? cuts::btag_2017::pt : cuts::btag_2016::pt;
}

double BTagger::EtaCut() const
{
    return period == analysis::Period::Run2017 ? cuts::btag_2017::eta : cuts::btag_2016::eta;
}

}
