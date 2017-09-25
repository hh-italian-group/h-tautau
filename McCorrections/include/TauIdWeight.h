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

namespace analysis {
namespace mc_corrections {


class TauIdWeight : public IWeightProvider {
public:
    using Event = ntuple::Event;
    using namespace boost::property_tree;

    struct Key{
        bool isData{false};
        bool isGenuineTau{false};
        analysis::DiscriminatorWP isowp;
        int decayMode;

        Key(bool _isData, bool _isGenuineTau, const std::string& _iso_wp, int _decayMode) :
            isData(_isData), isGenuineTau(_isGenuineTau), decayMode(_decayMode){
            isowp = Parse<analysis::DiscriminatorWP>(_iso_wp);
        }

        static Key Parse(const std::string& key_string)
        {
            Key key;
            std::map<std::string, bool> data_map = { {"data", true}, {"mc" , false } };
            std::map<std::string, bool> genuine_map = { {"genuine", true}, {"fake" , false } };

            std::vector<std::string> key_elements = analysis::SplitValueList(key_string,true, "_" );
            if(key_elements.size() != 4)
                throw exception("More than 4 key params found");
            for (unsigned n = 0; n < key_elements.size(); ++n){
                if(!data_map.count(key_elements.at(0)))
                    throw exception("It cannot be set if it is data or MC.");
                key.isData = data_map.at(key_elements.at(0));
                if(!genuine_map.count(key_elements.at(1)))
                    throw exception("It cannot be set if it is genuine or fake tau.");
                key.isGenuineTau = genuine_map.at(key_elements.at(1));
                std::string iso_string = key_elements.at(2).substr(0, key_elements.at(2) - 3);
                key.isowp = Parse<analysis::DiscriminatorWP>(iso_string);
                std::string idm_string = key_elements.at(3).substr(key_elements.at(2) - 2, 50);
                key.decayMode = Parse<int>(idm_string);
            }
            return key;
        }

    };

    struct Parameters{
        double alpha;
        double m_0;
        double sigma;
        double norm;
        double n;

        Parameters() :
            alpha(0.0), m_0(0.0), sigma(0.0), norm(0.0), n(0.0) {}

        Parameters(double _alpha, double _m_0, double _sigma, double _norm, double _n) :
            alpha(_alpha), m_0(_m_0), sigma(_sigma), norm(_norm), n(_n) {}

        static Parameters Parse(const ptree& pt)
        {
            Parameters parameters;
            for (auto& prop : pt.second){
                parameters.alpha = Parse<double>(prop.second.get_child("alpha"));
                parameters.m_0 = Parse<double>(prop.second.get_child("m_{0}"));
                parameters.sigma = Parse<double>(prop.second.get_child("sigma"));
                parameters.norm = Parse<double>(prop.second.get_child("norm"));
                parameters.n = Parse<double>(prop.second.get_child("n"));
            }
            return parameters;
        }
    };

    TauIdWeight(const std::string& tauId_input, const std::string& _iso_type) : iso_type(_iso_type)
    {
        ptree pt;
        json_parser::read_json(tauId_input, pt);

        for (auto& type : pt){
            std::cout << type.first << "\n";
            const Key key = Key::Parse(type.first);
            const Parameters parameters = Parameters::Parse(type.second);
            tauIdparam_map[key] = parameters;
        }

    }

    double EvaluateEfficiency(const double pt, const Parameters& parameters)
    {
        return analysis::crystalballEfficiency(pt,parameters.m_0,parameters.sigma,parameters.alpha,
                                               parameters.n,parameters.norm);
    }

    double GetWeight(const Event& event) const
    {
        const Channel channel = static_cast<Channel>(event.channelId);
        if(channel == Channel::TauTau) {
            Key key_data(true,true,"Medium",1);
            Key key_mc(false,true,"Medium",1);
            double eff_data_1, eff_data_2, eff_mc_1, eff_mc_2;
            auto iter_data = tauIdparam_map.find(key_data);
            if(iter_data != tauIdparam_map.end()){
                eff_data_1 = EvaluateEfficiency(event.p4_1.pt(),tauIdparam_map.at(key_data));
                eff_data_2 = EvaluateEfficiency(event.p4_2.pt(),tauIdparam_map.at(key_data));
            }
            throw exception("Key data not found");

            auto iter_mc = tauIdparam_map.find(key_mc);
            if(iter_mc != tauIdparam_map.end()){
                eff_mc_1 = EvaluateEfficiency(event.p4_1.pt(),tauIdparam_map.at(key_mc));
                eff_mc_2 = EvaluateEfficiency(event.p4_2.pt(),tauIdparam_map.at(key_mc));
            }
            throw exception("Key mc not found");
            double sf_1 = eff_data_1/eff_mc_1;
            double sf_2 = eff_data_2/eff_mc_2;
            return sf_1 * sf_2;
        }
        return 1.0;
    }

    virtual double Get(const Event& event) const override { return GetWeight(event); }

    virtual double Get(const ntuple::ExpressEvent& /*event*/) const override
    {
        throw exception("ExpressEvent is not supported in LeptonWeights::Get.");
    }

private:
    std::map<Key, Parameters> tauIdparam_map;
    std::string iso_type;

};

} // namespace mc_corrections
} // namespace analysis
