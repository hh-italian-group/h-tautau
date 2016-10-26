/*! MSSM Higgs to tautau selection.
Defined only cuts that are different from the H to tautau baseline selection.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once
namespace cuts {
  namespace H_tautau_2016_mssm {
    namespace MuTau {
      namespace muonID {
	constexpr double pt = 23;
      }
      namespace tauID { 
	constexpr double pt = 30; // >
      }
    }
    namespace ETau{
      namespace electronID {
	constexpr double pt = 26; // >
      } 
      namespace tauID { 
	constexpr double pt = 30; // >
      }
    }
  } // namespace H_tautau_mssm_2016
} // namespace cuts
