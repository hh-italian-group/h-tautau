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

using namespace boost::property_tree;

class TauIdWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;


    struct Key{
        bool isData{false};
        bool isGenuineTau{false};
        DiscriminatorWP isowp{DiscriminatorWP::Medium};
        int decayMode{0};

//        Key() {}
//        Key(bool _isData, bool _isGenuineTau, const std::string& _iso_wp, int _decayMode) :
//            isData(_isData), isGenuineTau(_isGenuineTau),isowp(Parse<analysis::DiscriminatorWP>(_iso_wp)),
//            decayMode(_decayMode) {}

        static Key Parse(const std::string& key_string)
        {
            Key key;
            static const std::map<std::string, bool> data_map = { {"data", true}, {"mc" , false } };
            static const std::map<std::string, bool> genuine_map = { {"genuine", true}, {"fake" , false } };

            std::vector<std::string> key_elements = SplitValueList(key_string,true, "_" );
            if(key_elements.size() != 4)
                throw exception("More than 4 key params found");

            if(!data_map.count(key_elements.at(0)))
                throw exception("Invalid tauId key %1%." ) %key_string;
            key.isData = data_map.at(key_elements.at(0));
            if(!genuine_map.count(key_elements.at(1)))
                throw exception("Invalid tauId key %1%." ) %key_string;
            key.isGenuineTau = genuine_map.at(key_elements.at(1));
            if(key_elements.at(2).substr(key_elements.at(2).size() - 3) != "Iso")
                throw exception("Invalid tauId key %1%." ) %key_string;
            if(!TryParse(key_elements.at(2).substr(0, key_elements.at(2).size() - 3), key.isowp))
                throw exception("Invalid tauId key %1%." ) %key_string;
            if(key_elements.at(3).substr(0,2) != "dm")
                throw exception("Invalid tauId key %1%." ) %key_string;
            if(!TryParse(key_elements.at(3).substr(2), key.decayMode))
                throw exception("Invalid tauId key %1%." ) %key_string;

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
        double alpha;
        double m_0;
        double sigma;
        double norm;
        double n;

//        Parameters() :
//            alpha(0.0), m_0(0.0), sigma(0.0), norm(0.0), n(0.0) {}

//        Parameters(double _alpha, double _m_0, double _sigma, double _norm, double _n) :
//            alpha(_alpha), m_0(_m_0), sigma(_sigma), norm(_norm), n(_n) {}

        static Parameters Parse(const ptree& pt)
        {
            Parameters parameters;


            parameters.alpha = pt.get_child("alpha").get_value<double>();
            parameters.m_0 = pt.get_child("m_{0}").get_value<double>();
            parameters.sigma = pt.get_child("sigma").get_value<double>();
            parameters.norm = pt.get_child("norm").get_value<double>();
            parameters.n = pt.get_child("n").get_value<double>();

            return parameters;
        }
    };

    TauIdWeight(const std::string& tauId_input, DiscriminatorWP _iso_wp) : iso_wp(_iso_wp)
    {
        ptree pt;
        json_parser::read_json(tauId_input, pt);

        for (auto& type : pt){
            std::cout << "Type first: " << type.first << "\n";
            const Key key = Key::Parse(type.first);
            std::cout << "Type second: " << type.second.get_value<std::string>() << "\n";
            tauIdparam_map[key] = Parameters::Parse(type.second);
        }

    }

    virtual double Get(const Event& event) const override
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel != Channel::TauTau) return 1.;

        const Key key_data{true, true, iso_wp, 1};
        const Key key_mc{false, true, iso_wp, 1};
        double eff_data_1, eff_data_2, eff_mc_1, eff_mc_2;
        auto iter_data = tauIdparam_map.find(key_data);
        if(iter_data == tauIdparam_map.end())
            throw exception("Key data not found for TauId weight.");
        eff_data_1 = EvaluateEfficiency(event.p4_1.pt(),tauIdparam_map.at(key_data));
        eff_data_2 = EvaluateEfficiency(event.p4_2.pt(),tauIdparam_map.at(key_data));

        auto iter_mc = tauIdparam_map.find(key_mc);
        if(iter_mc == tauIdparam_map.end())
            throw exception("Key mc not found for TauId weight.");
        eff_mc_1 = EvaluateEfficiency(event.p4_1.pt(),tauIdparam_map.at(key_mc));
        eff_mc_2 = EvaluateEfficiency(event.p4_2.pt(),tauIdparam_map.at(key_mc));

        double sf_1 = eff_data_1/eff_mc_1;
        double sf_2 = eff_data_2/eff_mc_2;
        return sf_1 * sf_2;

    }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in TauIdWeight::Get.");
    }

private:

    double EvaluateEfficiency(double pt, const Parameters& parameters) const
    {
        return analysis::crystalballEfficiency(pt,parameters.m_0,parameters.sigma,parameters.alpha,
                                               parameters.n,parameters.norm);
    }


    std::map<Key, Parameters> tauIdparam_map;
    DiscriminatorWP iso_wp;

};

} // namespace mc_corrections
} // namespace analysis
