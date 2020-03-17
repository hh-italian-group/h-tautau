/*! b-jet selection recommended by the b-tag POG.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

namespace cuts {
    namespace btag_Run2 {
        constexpr double pt = 20; // >
        constexpr double eta = 2.4; // <

        /*
        CSV definitions are taken from the b-tag POG 80X recommendation TWiki:
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco.
        DeepCSV and DeepFlavour definitions are taken from the b-tag 2016 Legacy recommendation Twiki:
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
        */
        namespace Run2016 {
            constexpr double CSVv2L = 0.5426; // >
            constexpr double CSVv2M = 0.8484; // >
            constexpr double CSVv2T = 0.9535; // >

            constexpr double DeepCSVL = 0.2217; // >
            constexpr double DeepCSVM = 0.6321; // >
            constexpr double DeepCSVT = 0.8953; // >

            constexpr double deepFlavourL = 0.0614; // >
            constexpr double deepFlavourM = 0.3093; // >
            constexpr double deepFlavourT = 0.7221; // >
        }

        /*
        All definitions are taken from the b-tag POG 94X recommendation TWiki:
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X.
        */
        namespace Run2017 {
            constexpr double CSVv2L = 0.5803; // >
            constexpr double CSVv2M = 0.8838; // >
            constexpr double CSVv2T = 0.9693; // >

            constexpr double deepCSVv2L = 0.1522; // >
            constexpr double deepCSVv2M = 0.4941; // >
            constexpr double deepCSVv2T = 0.8001; // >

            constexpr double deepFlavourL = 0.0521; // >
            constexpr double deepFlavourM = 0.3033; // >
            constexpr double deepFlavourT = 0.7489; // >
        }

        /*
        all definitions are taken from the b-tag POG 102X recommendation TWiki:
        https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X.
        */
        namespace Run2018 {
            constexpr double DeepCSVL = 0.1241; // >
            constexpr double DeepCSVM = 0.4184; // >
            constexpr double DeepCSVT = 0.7527; // >

            constexpr double deepFlavourL = 0.0494; // >
            constexpr double deepFlavourM = 0.2770; // >
            constexpr double deepFlavourT = 0.7264; // >
        }
    }
}
