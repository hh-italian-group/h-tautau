/*! Definiton of analysis::FlatEventInfo class.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <TLorentzVector.h>
#include <TMatrixD.h>

#include "FlatTree.h"
#include "AnalysisTypes.h"
#include "AnalysisTools/Core/include/exception.h"
#include "KinFit.h"
#include "Htautau_Summer13.h"


namespace analysis {

struct FlatEventInfo {
    typedef std::pair<size_t, size_t> BjetPair;

    static size_t NumberOfCombinationPairs(size_t n_bjets)
    {
        return n_bjets * (n_bjets - 1) / 2;
    }

    static size_t CombinationPairToIndex(const BjetPair& pair, size_t n_bjets)
    {
        const size_t min = std::min(pair.first, pair.second);
        const size_t max = std::max(pair.first, pair.second);
        if(n_bjets < 2 || min == max || min >= n_bjets || max >= n_bjets)
            throw exception("bad combination pair (") << pair.first << ", " << pair.second
                                                      << ") for n b-jets = " << n_bjets << ".";
        return max - 1 + min * (2 * n_bjets - 3 - min) / 2;
    }

    static BjetPair CombinationIndexToPair(size_t index, size_t n_bjets)
    {
        if(n_bjets < 2 || index >= n_bjets * (n_bjets - 1) / 2)
            throw exception("bad combination index = ") << index << " for n b-jets = " << n_bjets << ".";

        for(size_t min = 0;; ++min) {
            const size_t l = CombinationPairToIndex(BjetPair(min, n_bjets - 1), n_bjets);
            if(l >= index) {
                const size_t max = index + n_bjets - 1 - l;
                return BjetPair(min, max);
            }
        }
    }

    const ntuple::Flat* event;
    ntuple::EventType eventType;
    analysis::Channel channel;
    std::vector<TLorentzVector> lepton_momentums;
    std::vector<TLorentzVector> bjet_momentums;
    BjetPair selected_bjets;
    bool has_bjet_pair;
    TLorentzVector MET, Htt, Htt_MET, Hbb, resonance;
    TMatrixD MET_covariance;
    bool recalculate_mass_KinFit;
    analysis::kinematic_fit::four_body::FitResults fitResults;
    analysis::EventEnergyScale eventEnergyScale;
    double mva_BDT, mva_BDTD, mva_BDTMitFisher;


    FlatEventInfo(const ntuple::Flat& _event, const BjetPair& _selected_bjets, bool _recalculate_mass_KinFit)
        : event(&_event), eventType(static_cast<ntuple::EventType>(_event.eventType)),
          channel(static_cast<analysis::Channel>(_event.channel)),
          lepton_momentums(2), bjet_momentums(_event.pt_Bjets.size()), selected_bjets(_selected_bjets),
          has_bjet_pair(false), MET_covariance(2, 2), recalculate_mass_KinFit(_recalculate_mass_KinFit),
          eventEnergyScale(static_cast<analysis::EventEnergyScale>(_event.eventEnergyScale)),
          mva_BDT(-1), mva_BDTD(-1), mva_BDTMitFisher(-1)
    {
        lepton_momentums.at(0).SetPtEtaPhiM(event->pt_1, event->eta_1, event->phi_1, event->m_1);
        lepton_momentums.at(1).SetPtEtaPhiM(event->pt_2, event->eta_2, event->phi_2, event->m_2);

        for(size_t n = 0; n < bjet_momentums.size(); ++n)
            bjet_momentums.at(n).SetPtEtaPhiE(event->pt_Bjets.at(n), event->eta_Bjets.at(n), event->phi_Bjets.at(n),
                                              event->energy_Bjets.at(n));

        MET.SetPtEtaPhiM(event->mvamet, 0, event->mvametphi, 0);
        MET_covariance(0, 0) = event->mvacov00;
        MET_covariance(1, 0) = event->mvacov10;
        MET_covariance(0, 1) = event->mvacov01;
        MET_covariance(1, 1) = event->mvacov11;

        Htt = lepton_momentums.at(0) + lepton_momentums.at(1);
        Htt_MET = Htt + MET;
        has_bjet_pair = selected_bjets.first < bjet_momentums.size() && selected_bjets.second < bjet_momentums.size();
        if(has_bjet_pair) {
            Hbb = bjet_momentums.at(selected_bjets.first) + bjet_momentums.at(selected_bjets.second);
            resonance = Htt_MET + Hbb;
            if (recalculate_mass_KinFit){
                using namespace analysis::kinematic_fit;
                const four_body::FitInput four_body_input(bjet_momentums.at(selected_bjets.first),
                                                          bjet_momentums.at(selected_bjets.second),
                                                          lepton_momentums.at(0), lepton_momentums.at(1),
                                                          MET, MET_covariance);

                fitResults = Fit(four_body_input);
                if (fitResults.convergence == 0){
                    std::cout << "kin fit has convergence = 0! event = " << event->evt << std::endl;
                }
            } else {
                fitResults.convergence = event->kinfit_bb_tt_convergence;
                fitResults.chi2 = event->kinfit_bb_tt_chi2;
                fitResults.pull_balance = event->kinfit_bb_tt_pull_balance;
                fitResults.has_valid_mass = event->kinfit_bb_tt_convergence > 0;
                fitResults.mass = event->kinfit_bb_tt_mass;
            }

        }
    }
};

typedef std::shared_ptr<FlatEventInfo> FlatEventInfoPtr;

} // namespace analysis
