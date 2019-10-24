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
            {JetOrdering::Pt, {
                {DiscriminatorWP::Medium, cuts::btag_2016::pt}
            }}

        }}
    };

    if(!working_points.count(period))
        throw exception("Period %1% is not supported.") % period;
    auto base_ordering = ordering;
    if(ordering == JetOrdering::HHJetTag) base_ordering = JetOrdering::DeepFlavour;
    if(!working_points.at(period).count(base_ordering))
            throw exception("Jet Ordering %1% is not supported.") % base_ordering;

    cut = &working_points.at(period).at(base_ordering);
}

double BTagger::BTag(const ntuple::Event& event, size_t jet_index,
    analysis::UncertaintySource unc_source,analysis::UncertaintyScale unc_scale, bool use_base_ordering) const
{
    double tag = 0;
    if(use_base_ordering){
        if(ordering==JetOrdering::Pt) tag = event.jets_p4.at(jet_index).Pt();
        else if(ordering==JetOrdering::DeepCSV) tag = event.jets_deepCsv_BvsAll.at(jet_index);
        else if (ordering==JetOrdering::CSV) tag = event.jets_csv.at(jet_index);
        else if(ordering==JetOrdering::DeepFlavour) tag = (event.jets_deepFlavour_b.at(jet_index) +
                event.jets_deepFlavour_bb.at(jet_index) + event.jets_deepFlavour_lepb.at(jet_index));
    }
    else if(ordering==JetOrdering::HHJetTag) {
        for(unsigned n = 0; n < event.jet_hh_score_value.size(); ++n){
            if(event.jet_hh_score_index.at(n) != jet_index) continue;
            if(event.jet_hh_score_unc_source.at(n) != static_cast<int>(unc_source)) continue;
            if(event.jet_hh_score_unc_scale.at(n) != static_cast<int>(unc_scale)) continue;
            return event.jet_hh_score_value.at(n);
        }
    }
    else
        throw exception("Jet Ordering %1% is not supported") % ordering;
    return tag;
}

double BTagger::BTag(const ntuple::TupleJet& jet, analysis::UncertaintySource unc_source,
    analysis::UncertaintyScale unc_scale, bool use_base_ordering) const
{
    double tag = 0;
    if(use_base_ordering){
        if(ordering==JetOrdering::Pt) tag = jet.p4().Pt();
        else if(ordering==JetOrdering::DeepCSV) tag = jet.deepcsv();
        else if(ordering==JetOrdering::CSV) tag = jet.csv();
        else if(ordering==JetOrdering::DeepFlavour) tag = jet.deepFlavour();
    }
    else if(ordering==JetOrdering::HHJetTag) return jet.hh_tag(unc_source,unc_scale);
    else
        throw exception("Jet Ordering %1% is not supported") % ordering;
    return tag;
}

bool BTagger::Pass(const ntuple::Event& event, size_t jet_index,
    analysis::UncertaintySource unc_source, analysis::UncertaintyScale unc_scale, DiscriminatorWP wp) const
{
    if(!cut->count(wp))
        throw exception("Working point %1% is not supported.") %wp;
    else return BTag(event, jet_index,unc_source,unc_scale,true) > cut->at(wp);
}

bool BTagger::Pass(const ntuple::TupleJet& jet, analysis::UncertaintySource unc_source,
    analysis::UncertaintyScale unc_scale, DiscriminatorWP wp) const
{
    if(!cut->count(wp))
        throw exception("Working point %1% is not supported.") %wp;
    else return BTag(jet,unc_source,unc_scale,true) > cut->at(wp);
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
