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

    using ptree = boost::property_tree::ptree;

    struct Key {
        bool isData{false};
        bool isGenuineTau{false};
        DiscriminatorWP isowp{DiscriminatorWP::Medium};
        int decayMode{0};

        static Key Parse(const std::string& key_string)
        {
            Key key;
            static const std::map<std::string, bool> data_map = { {"data", true}, {"mc" , false } };
            static const std::map<std::string, bool> genuine_map = { {"genuine", true}, {"fake" , false } };

            std::vector<std::string> key_elements = SplitValueList(key_string, true, "_");
            if(key_elements.size() != 4)
                throw exception("More than 4 key params found");

            if(!data_map.count(key_elements.at(0)))
                throw exception("Invalid tauId key %1%." ) % key_string;
            key.isData = data_map.at(key_elements.at(0));
            if(!genuine_map.count(key_elements.at(1)))
                throw exception("Invalid tauId key %1%." ) % key_string;
            key.isGenuineTau = genuine_map.at(key_elements.at(1));
            if(key_elements.at(2).substr(key_elements.at(2).size() - 3) != "Iso")
                throw exception("Invalid tauId key %1%." ) % key_string;
            if(!TryParse(key_elements.at(2).substr(0, key_elements.at(2).size() - 3), key.isowp))
                throw exception("Invalid tauId key %1%." ) % key_string;
            if(key_elements.at(3).substr(0,2) != "dm")
                throw exception("Invalid tauId key %1%." ) % key_string;
            if(!TryParse(key_elements.at(3).substr(2), key.decayMode))
                throw exception("Invalid tauId key %1%." ) % key_string;

            return key;
        }

        bool operator < (const Key& other) const
        {
            if(isData != other.isData) return isData < other.isData;
            if(isGenuineTau != other.isGenuineTau) return isGenuineTau < other.isGenuineTau;
            if(isowp != other.isowp) return isowp < other.isowp;
            return decayMode < other.decayMode;
        }
    };

    struct Parameters{
        double alpha, m_0, sigma, norm, n;

        double EvaluateEfficiency(double pt) const
        {
            return pt < 1000 ? crystalball(pt, m_0, sigma, alpha, n, norm) : 1;
        }

        static Parameters Parse(const ptree& entry)
        {
            Parameters parameters;

            parameters.alpha = entry.get_child("alpha").get_value<double>();
            parameters.m_0 = entry.get_child("m_{0}").get_value<double>();
            parameters.sigma = entry.get_child("sigma").get_value<double>();
            parameters.norm = entry.get_child("norm").get_value<double>();
            parameters.n = entry.get_child("n").get_value<double>();

            return parameters;
        }
    };

    TauIdWeight2016(const std::string& tauId_input)
    {
        ptree property_tree;
        boost::property_tree::json_parser::read_json(tauId_input, property_tree);
        for (auto& type_entry : property_tree) {
            const Key key = Key::Parse(type_entry.first);
            tauIdparam_map[key] = Parameters::Parse(type_entry.second);
        }
    }
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch gen_match, int decay_mode,
                              DiscriminatorWP /*anti_ele_wp*/, DiscriminatorWP /*anti_mu_wp*/, DiscriminatorWP iso_wp) const override
    {
        return EvaluateSF(p4.pt(), gen_match, decay_mode, iso_wp);
    }

    private:
    double EvaluateSF(double pt, GenMatch gen_match, int decay_mode, DiscriminatorWP iso_wp) const
    {
        const bool is_genuine = gen_match == GenMatch::Tau;
        const double eff_data = EvaluateEfficiency(pt, true, is_genuine, decay_mode, iso_wp);
        const double eff_mc = EvaluateEfficiency(pt, false, is_genuine, decay_mode, iso_wp);
        return eff_data / eff_mc;
    }

    double EvaluateEfficiency(double pt, bool is_data, bool is_genuine, int decay_mode, DiscriminatorWP iso_wp) const
    {
        const Key key{is_data, is_genuine, iso_wp, decay_mode};
        auto iter = tauIdparam_map.find(key);
        if(iter == tauIdparam_map.end()) {
            const std::string name = is_data ? "data" : "mc";
            throw exception("Key %1% not found for TauId weight.") % name;
        }
        return iter->second.EvaluateEfficiency(pt);
    }

    std::map<Key, Parameters> tauIdparam_map;
};

class TauIdWeight2017 : public TauIdWeight {
public:
    virtual double GetIdIsoSF(const LorentzVectorM_Float& p4, GenMatch /*gen_match*/, int /*decay_mode*/, DiscriminatorWP anti_ele_wp,
                             DiscriminatorWP anti_mu_wp, DiscriminatorWP iso_wp) const override
    {
        auto tauSF = getTauIso(iso_wp);
        auto muonSF = getMuonMissIdSF(p4, anti_mu_wp);
        auto eleSF = getEleMissIdSF(p4, anti_ele_wp);

        return tauSF * muonSF * eleSF;
    }

private:
    double getMuonMissIdSF(const LorentzVectorM_Float& p4, DiscriminatorWP iso_wp) const {
        // https://indico.cern.ch/event/738043/contributions/3048471/attachments/1674773/2691664/TauId_26062018.pdf
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2017

        if(std::abs(p4.eta()) < 0.4){
            if(iso_wp == DiscriminatorWP::Loose)
                return 1.06;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 1.17;
            else throw exception("WP %1% is not supported.") % iso_wp;
        }
        else if(std::abs(p4.eta()) < 0.8 && std::abs(p4.eta()) > 0.4){
            if(iso_wp == DiscriminatorWP::Loose)
                return 1.02;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 1.29;
            else throw exception("WP %1% is not supported.") % iso_wp;
        }
        else if(std::abs(p4.eta()) > 0.8 && std::abs(p4.eta()) < 1.2){
            if(iso_wp == DiscriminatorWP::Loose)
                return 1.10;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 1.14;
            else throw exception("WP %1% is not supported.") % iso_wp;
        }
        else if(std::abs(p4.eta()) > 1.2 && std::abs(p4.eta()) < 1.7){
            if(iso_wp == DiscriminatorWP::Loose)
                return 1.03;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 0.93;
            else throw exception("WP %1% is not supported.") % iso_wp;
        }
        else if(std::abs(p4.eta()) < 1.7 && std::abs(p4.eta()) > 2.3){
            if(iso_wp == DiscriminatorWP::Loose)
                return 1.94;
            else if(iso_wp == DiscriminatorWP::Tight)
                return 1.67;
            else throw exception("WP %1% is not supported.") % iso_wp;
        }
        return 1;
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
        else
            return 1;

        return 1;
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

        return 1;
    }
};
} // namespace mc_corrections
} // namespace analysis
