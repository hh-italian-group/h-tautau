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

class TauIdWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;
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

    TauIdWeight(const std::string& tauId_input, DiscriminatorWP _iso_wp) : iso_wp(_iso_wp)
    {
        ptree property_tree;
        boost::property_tree::json_parser::read_json(tauId_input, property_tree);
        for (auto& type_entry : property_tree) {
            const Key key = Key::Parse(type_entry.first);
            tauIdparam_map[key] = Parameters::Parse(type_entry.second);
        }
    }

    virtual double Get(const Event& event) const override
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel != Channel::TauTau) return 1.;
        const double sf_1 = EvaluateSF(event.p4_1.pt(), static_cast<GenMatch>(event.gen_match_1), event.decayMode_1);
        const double sf_2 = EvaluateSF(event.p4_2.pt(), static_cast<GenMatch>(event.gen_match_2), event.decayMode_2);
        return sf_1 * sf_2;
    }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in TauIdWeight::Get.");
    }

private:
    double EvaluateSF(double pt, GenMatch gen_match, int decay_mode) const
    {
        const bool is_genuine = gen_match == GenMatch::Tau;
        const double eff_data = EvaluateEfficiency(pt, true, is_genuine, decay_mode);
        const double eff_mc = EvaluateEfficiency(pt, false, is_genuine, decay_mode);
        return eff_data / eff_mc;
    }

    double EvaluateEfficiency(double pt, bool is_data, bool is_genuine, int decay_mode) const
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
    DiscriminatorWP iso_wp;

};

} // namespace mc_corrections
} // namespace analysis
