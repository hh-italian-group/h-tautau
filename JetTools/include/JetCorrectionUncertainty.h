#pragma once

#include "SimpleJetCorrectionUncertainty.h"

namespace jec {

class JetCorrectionUncertainty
{
public:
    JetCorrectionUncertainty();
    JetCorrectionUncertainty(const std::string& fDataFile);
    JetCorrectionUncertainty(const JetCorrectorParameters& fParameters);
    ~JetCorrectionUncertainty();
    JetCorrectionUncertainty(const JetCorrectionUncertainty&) = delete;
    JetCorrectionUncertainty& operator= (const JetCorrectionUncertainty&) = delete;


    void setParameters(const std::string& fDataFile);
    void setAddLepToJet (bool fAddLepToJet);

    void setJetEta(float fEta);
    void setJetPt(float fPt);
    void setJetPhi(float fPhi);
    void setJetE(float fE);
    void setJetEMF(float fEMF);
    void setLepPx(float fPx);
    void setLepPy(float fPy);
    void setLepPz(float fPz);

    boost::optional<float> getUncertainty(bool fDirection);

private:

    //------------------------------------------------------------------------
    //--- Reads the parameter names and fills a vector of floats -------------
    //------------------------------------------------------------------------
    std::vector<float> fillVector(const std::vector<std::string>& fNames);

    //------------------------------------------------------------------------
    //--- Calculate the PtRel (needed for the SLB) ---------------------------
    //------------------------------------------------------------------------
    float getPtRel();

    //---- Member Data ---------
    float mJetE;
    float mJetEta;
    float mJetPt;
    float mJetPhi;
    float mJetEMF;
    float mLepPx;
    float mLepPy;
    float mLepPz;
    bool  mAddLepToJet;
    bool  mIsJetEset;
    bool  mIsJetPtset;
    bool  mIsJetPhiset;
    bool  mIsJetEtaset;
    bool  mIsJetEMFset;
    bool  mIsLepPxset;
    bool  mIsLepPyset;
    bool  mIsLepPzset;
    SimpleJetCorrectionUncertainty* mUncertainty;
};

}
