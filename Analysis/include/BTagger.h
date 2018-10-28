/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/TupleObjects.h"

namespace analysis {

enum class JetOrdering { NoOrdering, Pt, CSV, DeepCSV, DeepFlavour };
ENUM_NAMES(JetOrdering) = {
    { JetOrdering::NoOrdering, "NoOrdering" },
    { JetOrdering::Pt, "Pt" },
    { JetOrdering::CSV, "CSV" },
    { JetOrdering::DeepCSV, "DeepCSV" },
    { JetOrdering::DeepFlavour, "DeepFlavour" },
};

struct BTagger {
public:
    BTagger(Period period, JetOrdering _ordering) :
        ordering(_ordering)
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
        if(ordering != JetOrdering::Pt){
            if(!working_points.at(period).count(ordering))
                throw exception("Jet Ordering %1% is not supported.") % ordering;
        }

        cut = &working_points.at(period).at(ordering);
    }

    double BTag(const ntuple::Event& event, size_t jet_index) const
    {
        if(ordering==JetOrdering::Pt) return event.jets_p4.at(jet_index).Pt();
        else if(ordering==JetOrdering::DeepCSV) return event.jets_deepCsv_BvsAll.at(jet_index);
        else if (ordering==JetOrdering::CSV) return event.jets_csv.at(jet_index);
        else if(ordering==JetOrdering::DeepFlavour) return (event.jets_deepFlavour_b.at(jet_index) +
                event.jets_deepFlavour_bb.at(jet_index) + event.jets_deepFlavour_lepb.at(jet_index));
        else
            throw exception("Jet Ordering %1% is not supported") % ordering;
    }

    double BTag(const ntuple::TupleJet& jet) const
    {
        if(ordering==JetOrdering::Pt) return jet.p4().Pt();
        else if(ordering==JetOrdering::DeepCSV) return jet.deepcsv();
        else if(ordering==JetOrdering::CSV) return jet.csv();
        else if(ordering==JetOrdering::DeepFlavour) return jet.deepFlavour();
        else
            throw exception("Jet Ordering %1% is not supported") % ordering;
    }

    bool Pass(const ntuple::Event& event, size_t jet_index, DiscriminatorWP wp = DiscriminatorWP::Medium) const
    {
        if(!cut->count(wp))
            throw exception("Working point %1% is not supported.") %wp;
        else return BTag(event, jet_index) > cut->at(wp);
    }

    bool Pass(const ntuple::TupleJet& jet, DiscriminatorWP wp = DiscriminatorWP::Medium) const
    {
        if(!cut->count(wp))
            throw exception("Working point %1% is not supported.") %wp;
        else return BTag(jet) > cut->at(wp);
    }

private:
    JetOrdering ordering;
    const std::map<DiscriminatorWP,double>* cut;
};
  
}
