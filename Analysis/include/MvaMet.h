/*! Definition of wrapper for MVA MET code to apply MVA corrections for pfMet.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "METPUSubtraction/source/GBRForest.cxx"
#include "METPUSubtraction/source/GBRTree.cxx"
#include "METPUSubtraction/source/mvaMEtUtilities.cc"
#include "METPUSubtraction/source/PFMETAlgorithmMVA.cc"

#include "Candidate.h"
#include "Htautau_Summer13.h"

namespace analysis {

class MvaMetProducer {
public:
    MvaMetProducer(double dZcut, const std::string& inputFileNameU, const std::string& inputFileNameDPhi,
                   const std::string& inputFileNameCovU1, const std::string& inputFileNameCovU2)
        : metAlgo(dZcut), useType1Correction(false), minCorrJetPt(-1)
    {
        metAlgo.initialize(inputFileNameU, inputFileNameDPhi, inputFileNameCovU1, inputFileNameCovU2);
    }

    ntuple::MET ComputeMvaMet(const CandidatePtr& signalCandidate, const ntuple::PFCandVector& pfCandidates,
                              const ntuple::JetVector& jets, const VertexPtr& selectedVertex,
                              const VertexPtrVector& goodVertices)
    {
        const static bool debug = false;
        const auto leptonInfo = ComputeLeptonInfo(signalCandidate);
        auto pfCandidateInfo = ComputePFCandidateInfo(pfCandidates, selectedVertex->GetPosition());
        const auto vertexInfo = ComputeVertexInfo(goodVertices);
        const auto jetInfo = ComputeJetInfo(jets, leptonInfo, pfCandidateInfo);
        metAlgo.setInput(leptonInfo, jetInfo, pfCandidateInfo, vertexInfo);
        metAlgo.setHasPhotons(false);
        metAlgo.evaluateMVA();
        const TLorentzVector& metMomentum = metAlgo.getMEt();
        const TMatrix& metCov = metAlgo.getMEtCov();
        ntuple::MET mvaMET;
        mvaMET.phi = metMomentum.Phi();
        mvaMET.pt = metMomentum.Pt();
        mvaMET.significanceMatrix = ntuple::SignificanceMatrixToVector(metCov);
        if (debug){
            std::cout << "SIGNAL:\n";
            for ( const auto& lepton : leptonInfo ){
                std::cout << lepton.p4_ << ", charge fraction= " << lepton.chargedFrac_ << std::endl;
            }
            std::cout << "PFCandidates: " << pfCandidateInfo.size() << "\n";
            for ( const auto& pfCand : pfCandidateInfo ){
                std::cout << pfCand.p4_ << ", dZ= " << pfCand.dZ_ << std::endl;
            }
            std::cout << "Jets: " << jetInfo.size() << "\n";
            for ( const auto& jet : jetInfo ){
                std::cout << jet.p4_ << ", mva= " << jet.mva_ << ", neutralFrac= " << jet.neutralEnFrac_ << std::endl;
            }
            std::cout << "Vertices: " << vertexInfo.size() << "\n";
            for ( const auto& vertex : vertexInfo ){
                std::cout << vertex <<  std::endl;
            }
            std::cout << "MET:\n";
            std::cout << metMomentum << std::endl;
            std::cout << "Cov:\n";
            metAlgo.getMEtCov().Print();
            metAlgo.print(std::cout);
        }
        return mvaMET;
    }

    ntuple::MET ComputePFMet(const ntuple::PFCandVector& pfCandidates, const VertexPtr& selectedVertex)
    {
        auto pfCandidateInfo = ComputePFCandidateInfo(pfCandidates, selectedVertex->GetPosition());
        mvaMEtUtilities metUtilities;
        CommonMETData pfCandSum = metUtilities.computePFCandSum(pfCandidateInfo, 0.1, 2);
        const TVector2 vectorialMET(-pfCandSum.mex,-pfCandSum.mey);
        ntuple::MET pfMET;
        pfMET.sumEt = pfCandSum.sumet;
        pfMET.phi = vectorialMET.Phi();
        pfMET.pt = std::sqrt(vectorialMET.Px()*vectorialMET.Px() + vectorialMET.Py()*vectorialMET.Py());
        return pfMET;
    }

private:
    static double DefaultDeltaZ() { return -999.; }

    std::vector<mvaMEtUtilities::leptonInfo> ComputeLeptonInfo(const CandidatePtr& signalCandidate)
    {
        std::vector<mvaMEtUtilities::leptonInfo> leptonInfos;
        for(const auto& daughter : signalCandidate->GetFinalStateDaughters()) {
            mvaMEtUtilities::leptonInfo info;
            info.p4_ = daughter->GetMomentum();
            info.chargedFrac_ = ComputeChargedFraction(daughter);
            leptonInfos.push_back(info);
        }
        return leptonInfos;
    }

    std::vector<mvaMEtUtilities::pfCandInfo> ComputePFCandidateInfo(const ntuple::PFCandVector& pfCandidates,
                                                                    const TVector3& selectedVertex)
    {
        std::vector<mvaMEtUtilities::pfCandInfo> candInfos;
        for(const ntuple::PFCand& candidate : pfCandidates) {
            mvaMEtUtilities::pfCandInfo info;
            info.p4_.SetPtEtaPhiM(candidate.pt, candidate.eta, candidate.phi, candidate.mass); //with energy is the same
            if(candidate.haveTrackInfo) {
                const TVector3 trkV(candidate.trk_vx, candidate.trk_vy, candidate.trk_vz);
                TVector3 trkP;
                trkP.SetPtEtaPhi(candidate.pt, candidate.eta, candidate.phi);
                info.dZ_ = std::abs(Calculate_dz(trkV, selectedVertex, trkP));
            }
            else
                info.dZ_ = DefaultDeltaZ();
            candInfos.push_back(info);
        }
        return candInfos;
    }

    std::vector<TVector3> ComputeVertexInfo(const VertexPtrVector& goodVertices)
    {
        std::vector<TVector3> vertexInfos;
        for(const VertexPtr& vertex : goodVertices)
            vertexInfos.push_back(vertex->GetPosition());
        return vertexInfos;
    }

    std::vector<mvaMEtUtilities::JetInfo> ComputeJetInfo(const ntuple::JetVector& jets,
                                                         const std::vector<mvaMEtUtilities::leptonInfo>& signalLeptons,
                                                         std::vector<mvaMEtUtilities::pfCandInfo>& pfCandidates)
    {
        const static bool debug = false;
        static const double MinDeltaRtoSignalObjects = 0.5;
        static const double MinJetPtForPFcandCreation = 10.0;
        static const double MaxJetEtaForNeutralEnFrac = 2.5;

        std::vector<mvaMEtUtilities::JetInfo> jetInfos;
        for(const ntuple::Jet& jet : jets) {
            //if(!jet.passLooseID) continue;
            if(!cuts::Htautau_Summer13::jetID::passPFLooseId(jet)) continue;
            mvaMEtUtilities::JetInfo jetInfo;
            jetInfo.p4_.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass);
            double lType1Corr = 0;
            if(useType1Correction) {
                double pCorr = jet.correction;
                lType1Corr = jet.pt - pCorr * jet.pt_raw;
                TLorentzVector pType1Corr;
                pType1Corr.SetPtEtaPhiM(lType1Corr, 0, jet.phi, 0);
                bool pOnLepton = false;
                for(const mvaMEtUtilities::leptonInfo& signal : signalLeptons) {
                    if(signal.p4_.DeltaR(jetInfo.p4_) < MinDeltaRtoSignalObjects) pOnLepton = true;
                }
                if(jet.pt > MinJetPtForPFcandCreation && !pOnLepton) {
                    mvaMEtUtilities::pfCandInfo candInfo;
                    candInfo.p4_ = pType1Corr;
                    candInfo.dZ_ = DefaultDeltaZ();
                    pfCandidates.push_back(candInfo);
                }
                //lType1Corr = pCorr*jet.pt_raw - jet.pt_raw;
                lType1Corr /= jet.pt;
            }
            if(jetInfo.p4_.Pt() <= minCorrJetPt) continue;
            jetInfo.mva_ = jet.puIdMVA_met;
            jetInfo.neutralEnFrac_ = jet.neutralEmEnergyFraction + jet.neutralHadronEnergyFraction;
            if(std::abs(jet.eta) > MaxJetEtaForNeutralEnFrac) jetInfo.neutralEnFrac_ = 1.0;
//            if(useType1Correction) jetInfo.neutralEnFrac_ -= lType1Corr*jetInfo.neutralEnFrac_;
            if(useType1Correction) jetInfo.neutralEnFrac_ += lType1Corr;
            if (debug)
                std::cout << jetInfo.p4_ << ", mva= " << jetInfo.mva_ << ", neutralFrac= " << jetInfo.neutralEnFrac_
                          << ", pt_raw= " << jet.pt_raw << std::endl;

            jetInfos.push_back(jetInfo);
        }
        return jetInfos;
    }

    double ComputeChargedFraction(const CandidatePtr& candidate)
    {
        if(candidate->GetType() == Candidate::Type::Muon || candidate->GetType() == Candidate::Type::Electron)
            return 1.0;
        if(candidate->GetType() != Candidate::Type::Tau)
            throw exception("Unsupported candidate type to compute charged fraction.");
        const ntuple::Tau& tau = candidate->GetNtupleObject<ntuple::Tau>();
        double ptTotal = 0.0, ptCharged = 0.0;
        for(double pt : tau.signalChHadCand_Pt) {
            ptCharged += pt;
            ptTotal += pt;
        }
        for(double pt : tau.signalNeutrHadCand_Pt)
            ptTotal += pt;
        for(double pt : tau.signalGammaCand_Pt)
            ptTotal += pt;
        if(!ptTotal) ptTotal = 1.0;
        return ptCharged / ptTotal;
    }

private:
    PFMETAlgorithmMVA metAlgo;
    bool useType1Correction;
    double minCorrJetPt;
};

} // analysis
