/*! b-jet tagging.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

// put includes here
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Cuts/include/Btag_2017.h"
#include "h-tautau/Cuts/include/Btag_2016.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "h-tautau/Analysis/include/TupleObjects.h"

namespace analysis {

enum class JetOrdering { NoOrdering, Pt, CSV, DeepCSV, DeepFalvour };
ENUM_NAMES(JetOrdering) = {
    { JetOrdering::NoOrdering, "NoOrdering" },
    { JetOrdering::Pt, "Pt" },
    { JetOrdering::CSV, "CSV" },
    { JetOrdering::DeepCSV, "DeepCSV" },
    { JetOrdering::DeepFalvour, "DeepFalvour" },
};

struct BTagger {
public:
    BTagger(Period period, JetOrdering _ordering, DiscriminatorWP wp) :
        ordering(_ordering)
    {
        static const std::map<Period, std::map<JetOrdering, std::map<DiscriminatorWP, double>>> working_points = {
            {Period::Run2017, {
                {JetOrdering::DeepCSV, {
                    {DiscriminatorWP::Loose, cuts::btag_2017::deepCSVv2L},{DiscriminatorWP::Medium, cuts::btag_2017::deepCSVv2M},
                    {DiscriminatorWP::Tight, cuts::btag_2017::deepCSVv2T}}},
                {JetOrdering::CSV, {
                    {DiscriminatorWP::Loose, cuts::btag_2017::CSVv2L}, {DiscriminatorWP::Medium, cuts::btag_2017::CSVv2M},
                    {DiscriminatorWP::Tight, cuts::btag_2017::CSVv2T}}},
                {JetOrdering::DeepFalvour, {
                    {DiscriminatorWP::Loose, cuts::btag_2017::deepFlavourL},
                    {DiscriminatorWP::Medium, cuts::btag_2017::deepFlavourM},
                    {DiscriminatorWP::Tight, cuts::btag_2017::deepFlavourT}}}}},
            {Period::Run2016,
                {JetOrdering::CSV, {
                    {DiscriminatorWP::Loose, cuts::btag_2016::CSVv2L},
                    {DiscriminatorWP::Medium, cuts::btag_2016::CSVv2M},
                    {DiscriminatorWP::Tight, cuts::btag_2016::CSVv2T}}}}};

        if(!working_points.at(period).at(ordering).at(wp))
            throw exception("Working point %1% is not supported.") % wp;
        
        // need to do checks here that entry exists, otherwise throw an exception.
        cut = working_points.at(period).at(ordering).at(wp);
    }

    double BTag(const ntuple::Event& event, size_t jet_index) const
    {
        // ...
        if(ordering==JetOrdering::DeepCSV) return event.jets_deepCsv_BvsAll.at(jet_index);
        else if (ordering==JetOrdering::CSV) return event.jets_csv.at(jet_index);
        else if(ordering==JetOrdering::DeepFalvour) return (event.jets_deepFlavour_b.at(jet_index) +
                event.jets_deepFlavour_bb.at(jet_index) + event.jets_deepFlavour_lepb.at(jet_index));
    }

    double BTag(const ntuple::TupleJet& jet) const
    {
        if(ordering==JetOrdering::DeepCSV) return jet.deepcsv();
        else if(ordering==JetOrdering::CSV) return jet.csv();
        else if(ordering==JetOrdering::DeepFalvour) return jet.deepFlavour();
    }

    bool Pass(const ntuple::Event& event, size_t jet_index) const
    {
        return BTag(event, jet_index) > cut;
    }

    bool Pass(const ntuple::TupleJet& jet) const
    {
        return BTag(jet) > cut;
    }

private:
    JetOrdering ordering;
    double cut;
};
  
}
