/*! Various lepton weights.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "AnalysisTools/Core/include/AnalysisMath.h"
#include "h-tautau/Analysis/include/AnalysisTypes.h"
#include "h-tautau/Analysis/include/EventTuple.h"
#include "AnalysisTools/Core/include/RootExt.h"
#include "WeightProvider.h"
#include "TextIO.h"
#include "TString.h"

namespace analysis {
namespace mc_corrections {
//namespace detail{

class TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch gen_match, int decay_mode, DiscriminatorWP anti_ele_wp,
                      DiscriminatorWP anti_mu_wp, DiscriminatorWP iso_wp) const = 0;

    virtual double GetTauIdEfficiencyUncertainty(DiscriminatorWP iso_wp) const = 0;

    virtual double GetMuonMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP anti_mu_wp) const = 0;

    virtual double GetEleMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP anti_ele_wp) const = 0;

    virtual ~TauIdWeight() {}
};

class TauIdWeight2016 : public TauIdWeight {
public:

    virtual double GetIdIsoSF(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/, int /*decay_mode*/, DiscriminatorWP /*anti_ele_wp*/,
                             DiscriminatorWP /*anti_mu_wp*/, DiscriminatorWP /*iso_wp*/) const override
    {
        return 1;
    }

    virtual double GetTauIdEfficiencyUncertainty(DiscriminatorWP /*iso_wp*/) const override
    {
        return 0;
    }

    virtual double GetMuonMissIdUncertainty(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/, DiscriminatorWP /*anti_mu_wp*/) const override
    {
        return 0;
    }
    virtual double GetEleMissIdUncertainty(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/, DiscriminatorWP /*anti_ele_wp*/) const override
    {
        return 0;
    }


};

class TauIdWeight2017 : public TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch gen_match, int /*decay_mode*/, DiscriminatorWP anti_ele_wp,
                             DiscriminatorWP anti_mu_wp, DiscriminatorWP iso_wp) const override
    {
        auto tauSF = gen_match == GenMatch::Tau ? getTauIso(iso_wp).GetValue() : 1;
        auto muonSF = getMuonMissId(p4, gen_match, anti_mu_wp).GetValue();
        auto eleSF = getEleMissId(p4, gen_match, anti_ele_wp).GetValue();

        return tauSF * muonSF * eleSF;
    }

    virtual double GetTauIdEfficiencyUncertainty(DiscriminatorWP iso_wp) const override
    {
        return getTauIso(iso_wp).GetRelativeStatisticalError();
    }

    virtual double GetMuonMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP anti_mu_wp) const override
    {
        return getMuonMissId(p4, gen_match, anti_mu_wp).GetRelativeStatisticalError();
    }

    virtual double  GetEleMissIdUncertainty(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP anti_ele_wp) const override
    {
        return getEleMissId(p4, gen_match, anti_ele_wp).GetRelativeStatisticalError();
    }

private:

    PhysicalValue getMuonMissId(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP iso_wp) const {
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV

        if (!(gen_match == GenMatch::Muon || gen_match == GenMatch::TauMuon))
            return PhysicalValue(1,0);

        if(!(iso_wp == DiscriminatorWP::Loose || iso_wp == DiscriminatorWP::Tight))
            throw exception("WP %1% is not supported.") % iso_wp;

        if(std::abs(p4.eta()) < 0.4)
            return iso_wp == DiscriminatorWP::Loose ? PhysicalValue(1.06,0.05) : PhysicalValue(1.17,0.12);

        else if(std::abs(p4.eta()) >= 0.4 && std::abs(p4.eta()) < 0.8)
            return iso_wp == DiscriminatorWP::Loose ? PhysicalValue(1.02,0.04) : PhysicalValue(1.29,0.30);

        else if(std::abs(p4.eta()) >= 0.8 && std::abs(p4.eta()) < 1.2)
            return iso_wp == DiscriminatorWP::Loose ? PhysicalValue(1.10,0.04) : PhysicalValue(1.14,0.05);

        else if(std::abs(p4.eta()) >= 1.2 && std::abs(p4.eta()) < 1.7)
            return iso_wp == DiscriminatorWP::Loose ? PhysicalValue(1.03,0.18) : PhysicalValue(0.93,0.60);

        else if(std::abs(p4.eta()) >=1.7 && std::abs(p4.eta()) < 2.3)
                    return iso_wp == DiscriminatorWP::Loose ? PhysicalValue(1.94,0.35) : PhysicalValue(1.61,0.60);

        else throw exception("eta out of range");
    }

    PhysicalValue getEleMissId(const LorentzVectorM_Float& p4, GenMatch gen_match, DiscriminatorWP iso_wp) const{
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV

        if (!(gen_match == GenMatch::Electron || gen_match == GenMatch::TauElectron))
            return PhysicalValue(1,0);

        //Barrel ( abs(eta) < 1.460)
        if(std::abs(p4.eta()) < 1.460){
            if(iso_wp == DiscriminatorWP::VLoose)
                return PhysicalValue(1.09,0.01);
            else if(iso_wp == DiscriminatorWP::Loose)
                return PhysicalValue(1.17,0.04);
            else if(iso_wp == DiscriminatorWP::Medium)
                return  PhysicalValue(1.40,0.12);
            else if(iso_wp == DiscriminatorWP::Tight)
                return  PhysicalValue(1.80,0.20);
            else if(iso_wp == DiscriminatorWP::VTight)
                return PhysicalValue(1.96,0.27);

            else throw exception("WP %1% is not supported.") % iso_wp;
        }

        // Endcaps ( abs(eta) > 1.558)
        else if(std::abs(p4.eta()) > 1.558){
            if(iso_wp == DiscriminatorWP::VLoose)
                return PhysicalValue(1.19,0.01);
            else if(iso_wp == DiscriminatorWP::Loose)
                return PhysicalValue(1.25,0.06);
            else if(iso_wp == DiscriminatorWP::Medium)
                return PhysicalValue(1.21,0.26);
            else if(iso_wp == DiscriminatorWP::Tight)
                return PhysicalValue(1.53,0.60);
            else if(iso_wp == DiscriminatorWP::VTight)
                return PhysicalValue(1.66,0.80);

            else throw exception("WP %1% is not supported.") % iso_wp;
        }

        //Gap between barrel and endcaps
        else if(std::abs(p4.eta()) >= 1.460 && std::abs(p4.eta()) <= 1.558)
            return PhysicalValue(1,0);

        else throw exception("eta out of range");
    }

    //Isolation sum with deltaR=0.5
    PhysicalValue getTauIso(DiscriminatorWP iso_wp) const{
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf

        if(iso_wp == DiscriminatorWP::VLoose)
            return  PhysicalValue(0.88,0.03);
        else if(iso_wp == DiscriminatorWP::Loose || iso_wp == DiscriminatorWP::Medium || iso_wp == DiscriminatorWP::Tight)
            return  PhysicalValue(0.89,0.03);
        else if(iso_wp == DiscriminatorWP::VTight)
            return  PhysicalValue(0.86,0.03);
        else if(iso_wp == DiscriminatorWP::VVTight)
            return  PhysicalValue(0.84,0.03);
        else throw exception("WP %1% is not supported.") % iso_wp;
    }
};
} // namespace mc_corrections
} // namespace analysis
