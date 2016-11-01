/*! SM Higgs to tautau selection.
Defined only cuts that are different from the H to tautau baseline selection.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
namespace cuts {
  namespace H_tautau_2016_sm {
    namespace MuTau {
      namespace muonID {
	constexpr double pt = 23;
        constexpr double eta = 2.4; // <
      }
    }
    namespace ETau{
      namespace electronID {
	constexpr double pt = 26; // >
      }
    }
  }// namespace H_tau_tau_2016_sm
}//namespace cuts 
