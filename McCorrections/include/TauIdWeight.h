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

    virtual ~TauIdWeight() {}
};

class TauIdWeight2016 : public TauIdWeight {
public:

    virtual double GetIdIsoSF(const LorentzVectorM_Float& /*p4*/, GenMatch /*gen_match*/, int /*decay_mode*/, DiscriminatorWP /*anti_ele_wp*/,
                             DiscriminatorWP /*anti_mu_wp*/, DiscriminatorWP /*iso_wp*/) const override
    {
        return 1;
    }
};

class TauIdWeight2017 : public TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch gen_match, int /*decay_mode*/, DiscriminatorWP anti_ele_wp,
                             DiscriminatorWP anti_mu_wp, DiscriminatorWP iso_wp) const override
    {
        auto tauSF = getTauIso(iso_wp);
        auto muonSF = (gen_match == GenMatch::Muon || gen_match == GenMatch::TauMuon) ? getMuonMissIdSF(p4, anti_mu_wp) : 1;
        auto eleSF = (gen_match == GenMatch::Electron || gen_match == GenMatch::TauElectron) ? getEleMissIdSF(p4, anti_ele_wp) : 1;

        return tauSF * muonSF * eleSF;
    }

    virtual double GetIdIsoSFUncertainties(const LorentzVectorM_Float& p4, GenMatch gen_match, int /*decay_mode*/, DiscriminatorWP anti_ele_wp,
                             DiscriminatorWP anti_mu_wp, DiscriminatorWP iso_wp) const override
    {
        auto tauSFUnc = getTauIsoUncertainties(iso_wp);
        auto muonSFUnc = (gen_match == GenMatch::Muon || gen_match == GenMatch::TauMuon) ? getMuonMissIdSFUncertainties(p4, anti_mu_wp) : 0;
        auto eleSFUnc = (gen_match == GenMatch::Electron || gen_match == GenMatch::TauElectron) ? getEleMissIdSF(p4, anti_ele_wp) : 0;

        return std::sqrt(tauSFUnc*tauSFUnc+muonSFUnc*muonSFUnc+eleSFUnc*eleSFUnc);
    }

private:
    double getMuonMissIdSF(const LorentzVectorM_Float& p4, DiscriminatorWP iso_wp) const {
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017

        if(!(iso_wp == DiscriminatorWP::Loose || iso_wp == DiscriminatorWP::Tight))
            throw exception("WP %1% is not supported.") % iso_wp;

        if(std::abs(p4.eta()) < 0.4)
            return iso_wp == DiscriminatorWP::Loose ? 1.06 : 1.17;

        else if(std::abs(p4.eta()) >= 0.4 && std::abs(p4.eta()) < 0.8)
            return iso_wp == DiscriminatorWP::Loose ? 1.02 : 1.29;

        else if(std::abs(p4.eta()) >= 0.8 && std::abs(p4.eta()) < 1.2)
            return iso_wp == DiscriminatorWP::Loose ? 1.10 : 1.14;

        else if(std::abs(p4.eta()) >= 1.2 && std::abs(p4.eta()) < 1.7)
            return iso_wp == DiscriminatorWP::Loose ? 1.03 : 0.93;

        else if(std::abs(p4.eta()) >=1.7 && std::abs(p4.eta()) < 2.3)
                    return iso_wp == DiscriminatorWP::Loose ? 1.94 : 1.61;

        else throw exception("eta out of range");
    }

    double getEleMissIdSF(const LorentzVectorM_Float& p4, DiscriminatorWP iso_wp) const{
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        //Recommendation for Ele SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations

        //Barrel ( abs(eta) < 1.460)
        if(std::abs(p4.eta()) < 1.460){
            if(iso_wp == DiscriminatorWP::VLoose)
                return 1.09;
            else if(iso_wp == DiscriminatorWP::Loose)
                return 1.17;
            else if(iso_wp == DiscriminatorWP::Medium)
                return  1.40;
            else if(iso_wp == DiscriminatorWP::Tight)
                return  1.80;
            else if(iso_wp == DiscriminatorWP::VTight)
                return 1.96;

            else throw exception("WP %1% is not supported.") % iso_wp;
        }

        // Endcaps ( abs(eta) > 1.558)
        else if(std::abs(p4.eta()) > 1.558){
            if(iso_wp == DiscriminatorWP::VLoose)
                return 1.19;
            else if(iso_wp == DiscriminatorWP::Loose)
                return 1.25;
            else if(iso_wp == DiscriminatorWP::Medium)
                return 1.21;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 1.53;
            else if(iso_wp == DiscriminatorWP::VTight)
                return 1.66;

            else throw exception("WP %1% is not supported.") % iso_wp;
        }

        //Gap between barrel and endcaps
        else if(std::abs(p4.eta()) >= 1.460 && std::abs(p4.eta()) <= 1.558)
            return 1;

        else throw exception("eta out of range");
    }

    //Isolation sum with deltaR=0.5
    double getTauIso(DiscriminatorWP iso_wp) const{
        //Recommendation for Tau SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        // double tauIsoSF = 1;
        if(iso_wp == DiscriminatorWP::Medium)
            return  0.94;
        else if (iso_wp == DiscriminatorWP::Tight)
            return  0.93;
        else throw exception("WP %1% is not supported.") % iso_wp;
    }


    double getMuonMissIdSFUncertainties(const LorentzVectorM_Float& p4, DiscriminatorWP iso_wp) const {
            //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
            //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017

            if(!(iso_wp == DiscriminatorWP::Loose || iso_wp == DiscriminatorWP::Tight))
                throw exception("WP %1% is not supported.") % iso_wp;

            if(std::abs(p4.eta()) < 0.4)
                return iso_wp == DiscriminatorWP::Loose ? 0.05 : 0.12;

            else if(std::abs(p4.eta()) >= 0.4 && std::abs(p4.eta()) < 0.8)
                return iso_wp == DiscriminatorWP::Loose ? 0.04 : 0.30;

            else if(std::abs(p4.eta()) >= 0.8 && std::abs(p4.eta()) < 1.2)
                return iso_wp == DiscriminatorWP::Loose ? 0.04 : 0.05;

            else if(std::abs(p4.eta()) >= 1.2 && std::abs(p4.eta()) < 1.7)
                return iso_wp == DiscriminatorWP::Loose ? 1.03 : 0.93;

            else if(std::abs(p4.eta()) >=1.7 && std::abs(p4.eta()) < 2.3)
                        return iso_wp == DiscriminatorWP::Loose ? 1.94 : 1.61;

            else throw exception("eta out of range");
        }

    double getEleMissIdSFUncertainties(const LorentzVectorM_Float& p4, DiscriminatorWP iso_wp) const{
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        //Recommendation for Ele SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/Egamma2017DataRecommendations

        //Barrel ( abs(eta) < 1.460)
        if(std::abs(p4.eta()) < 1.460){
            if(iso_wp == DiscriminatorWP::VLoose)
                return 0.01;
            else if(iso_wp == DiscriminatorWP::Loose)
                return 0.04;
            else if(iso_wp == DiscriminatorWP::Medium)
                return  0.12;
            else if(iso_wp == DiscriminatorWP::Tight)
                return  0.20;
            else if(iso_wp == DiscriminatorWP::VTight)
                return 0.27;

            else throw exception("WP %1% is not supported.") % iso_wp;
        }

        // Endcaps ( abs(eta) > 1.558)
        else if(std::abs(p4.eta()) > 1.558){
            if(iso_wp == DiscriminatorWP::VLoose)
                return 0.01;
            else if(iso_wp == DiscriminatorWP::Loose)
                return 0.06;
            else if(iso_wp == DiscriminatorWP::Medium)
                return 0.26;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 0.60;
            else if(iso_wp == DiscriminatorWP::VTight)
                return 0.80;

            else throw exception("WP %1% is not supported.") % iso_wp;
        }

        //Gap between barrel and endcaps
        else if(std::abs(p4.eta()) >= 1.460 && std::abs(p4.eta()) <= 1.558)
            return 0;

        else throw exception("eta out of range");
    }

    //Isolation sum with deltaR=0.5
    double getTauIsoUncertainties(DiscriminatorWP iso_wp) const{
        //Recommendation for Tau SF https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
        //https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        if(iso_wp == DiscriminatorWP::Medium || iso_wp == DiscriminatorWP::Tight)
            return  0.05;
        else throw exception("WP %1% is not supported.") % iso_wp;
    }

};
} // namespace mc_corrections
} // namespace analysis
