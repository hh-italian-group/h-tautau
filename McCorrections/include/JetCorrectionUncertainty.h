#pragma once

#include "h-tautau/McCorrections/include/Serializable.h"
#include "h-tautau/McCorrections/include/JetCorrectorParameters.h"
#include "h-tautau/McCorrections/include/SimpleJetCorrectionUncertainty.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/Vector3D.h"
#include "Math/LorentzVector.h"
#include <vector>
#include <string>

namespace jec {

    class SimpleJetCorrectionUncertainty;
    class JetCorrectorParameters;

    class JetCorrectionUncertainty
    {
          public:
            JetCorrectionUncertainty() {
                mJetEta = -9999;
                mJetPt  = -9999;
                mJetPhi = -9999;
                mJetE   = -9999;
                mJetEMF = -9999;
                mLepPx  = -9999;
                mLepPy  = -9999;
                mLepPz  = -9999;
                mIsJetEset   = false;
                mIsJetPtset  = false;
                mIsJetPhiset = false;
                mIsJetEtaset = false;
                mIsJetEMFset = false;
                mIsLepPxset  = false;
                mIsLepPyset  = false;
                mIsLepPzset  = false;
                mAddLepToJet = false;
                mUncertainty = new SimpleJetCorrectionUncertainty();
            }

            JetCorrectionUncertainty(const std::string& fDataFile) {
                mJetEta = -9999;
                mJetPt  = -9999;
                mJetPhi = -9999;
                mJetE   = -9999;
                mJetEMF = -9999;
                mLepPx  = -9999;
                mLepPy  = -9999;
                mLepPz  = -9999;
                mIsJetEset   = false;
                mIsJetPtset  = false;
                mIsJetPhiset = false;
                mIsJetEtaset = false;
                mIsJetEMFset = false;
                mIsLepPxset  = false;
                mIsLepPyset  = false;
                mIsLepPzset  = false;
                mAddLepToJet = false;
                mUncertainty = new SimpleJetCorrectionUncertainty(fDataFile);
            }

            JetCorrectionUncertainty(const JetCorrectorParameters& fParameters) {
                mJetEta = -9999;
                mJetPt  = -9999;
                mJetPhi = -9999;
                mJetE   = -9999;
                mJetEMF = -9999;
                mLepPx  = -9999;
                mLepPy  = -9999;
                mLepPz  = -9999;
                mIsJetEset   = false;
                mIsJetPtset  = false;
                mIsJetPhiset = false;
                mIsJetEtaset = false;
                mIsJetEMFset = false;
                mIsLepPxset  = false;
                mIsLepPyset  = false;
                mIsLepPzset  = false;
                mAddLepToJet = false;
                mUncertainty = new SimpleJetCorrectionUncertainty(fParameters);
            }

            ~JetCorrectionUncertainty(){
                delete mUncertainty;
            }

            void setParameters  (const std::string& fDataFile){
                //---- delete the mParameters pointer before setting the new address ---
                delete mUncertainty;
                mUncertainty = new SimpleJetCorrectionUncertainty(fDataFile);
            }

            void setAddLepToJet (bool fAddLepToJet) {mAddLepToJet = fAddLepToJet;}

            //------------------------------------------------------------------------
            //--- Setters ------------------------------------------------------------
            //------------------------------------------------------------------------
            void setJetEta(float fEta)
            {
              mJetEta = fEta;
              mIsJetEtaset = true;
            }
            //------------------------------------------------------------------------
            void setJetPt(float fPt)
            {
              mJetPt = fPt;
              mIsJetPtset  = true;
            }
            //------------------------------------------------------------------------
            void setJetPhi(float fPhi)
            {
              mJetPhi = fPhi;
              mIsJetPhiset  = true;
            }
            //------------------------------------------------------------------------
            void setJetE(float fE)
            {
              mJetE = fE;
              mIsJetEset   = true;
            }
            //------------------------------------------------------------------------
            void setJetEMF(float fEMF)
            {
              mJetEMF = fEMF;
              mIsJetEMFset = true;
            }
            //------------------------------------------------------------------------
            void setLepPx(float fPx)
            {
              mLepPx = fPx;
              mIsLepPxset  = true;
            }
            //------------------------------------------------------------------------
            void setLepPy(float fPy)
            {
              mLepPy = fPy;
              mIsLepPyset  = true;
            }
            //------------------------------------------------------------------------
            void setLepPz(float fPz)
            {
              mLepPz = fPz;
              mIsLepPzset  = true;
            }
            //------------------------------------------------------------------------

            float getUncertainty(bool fDirection)
            {
                float result;
                std::vector<float> vx,vy;
                vx = fillVector(mUncertainty->parameters().definitions().binVar());
                vy = fillVector(mUncertainty->parameters().definitions().parVar());
                result = mUncertainty->uncertainty(vx,vy[0],fDirection);
                mIsJetEset   = false;
                mIsJetPtset  = false;
                mIsJetPhiset = false;
                mIsJetEtaset = false;
                mIsJetEMFset = false;
                mIsLepPxset  = false;
                mIsLepPyset  = false;
                mIsLepPzset  = false;
                return result;
            }

         private:
              JetCorrectionUncertainty(const JetCorrectionUncertainty&);
              JetCorrectionUncertainty& operator= (const JetCorrectionUncertainty&);

              //------------------------------------------------------------------------
              //--- Reads the parameter names and fills a vector of floats -------------
              //------------------------------------------------------------------------
              std::vector<float> fillVector(const std::vector<std::string>& fNames)
              {
                std::vector<float> result;
                for(unsigned i=0;i<fNames.size();i++)
                  {
                    if (fNames[i] == "JetEta")
                      {
                        if (!mIsJetEtaset) {
                      std::cout << "JetCorrectionUncertainty: jet eta is not set" << std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mJetEta);
                    }
                      }
                    else if (fNames[i] == "JetPt")
                      {
                        if (!mIsJetPtset){
                      std::cout << "JetCorrectionUncertainty: jet pt is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mJetPt);
                    }
                      }
                    else if (fNames[i] == "JetPhi")
                      {
                        if (!mIsJetPhiset) {
                      std::cout << "JetCorrectionUncertainty: jet phi is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mJetPt);
                    }
                      }
                    else if (fNames[i] == "JetE")
                      {
                        if (!mIsJetEset) {
                      std::cout << "JetCorrectionUncertainty: jet energy is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mJetE);
                    }
                      }
                    else if (fNames[i] == "JetEMF")
                      {
                        if (!mIsJetEMFset) {
                      std::cout << "JetCorrectionUncertainty: jet emf is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mJetEMF);
                    }
                      }
                    else if (fNames[i] == "LepPx")
                      {
                        if (!mIsLepPxset){
                      std::cout << "JetCorrectionUncertainty: lepton px is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mLepPx);
                    }
                      }
                    else if (fNames[i] == "LepPy")
                      {
                        if (!mIsLepPyset){
                      std::cout << "JetCorrectionUncertainty: lepton py is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mLepPy);
                    }
                      }
                    else if (fNames[i] == "LepPz")
                      {
                        if (!mIsLepPzset){
                      std::cout << "JetCorrectionUncertainty: lepton pz is not set"<< std::endl;
                      result.push_back(-999.0);
                    } else {
                      result.push_back(mLepPz);
                    }
                      }

                    else {
                  std::cout << "JetCorrectionUncertainty: unknown parameter "<<fNames[i] << std::endl;
                  result.push_back(-999.0);
                    }
                  }
                return result;
              }

              //------------------------------------------------------------------------
              //--- Calculate the PtRel (needed for the SLB) ---------------------------
              //------------------------------------------------------------------------
              float getPtRel()
              {
                typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<float> > PtEtaPhiELorentzVector;
                typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float> > XYZVector;
                PtEtaPhiELorentzVector jet;
                XYZVector lep;
                jet.SetPt(mJetPt);
                jet.SetEta(mJetEta);
                jet.SetPhi(mJetPhi);
                jet.SetE(mJetE);
                lep.SetXYZ(mLepPx,mLepPy,mLepPz);
                float lj_x = (mAddLepToJet) ? lep.X()+jet.Px() : jet.Px();
                float lj_y = (mAddLepToJet) ? lep.Y()+jet.Py() : jet.Py();
                float lj_z = (mAddLepToJet) ? lep.Z()+jet.Pz() : jet.Pz();
                // absolute values squared
                float lj2  = lj_x*lj_x+lj_y*lj_y+lj_z*lj_z;
                float pTrel2 = -999.0;
                if (lj2 > 0) {
                  float lep2 = lep.X()*lep.X()+lep.Y()*lep.Y()+lep.Z()*lep.Z();
                  // projection vec(mu) to lepjet axis
                  float lepXlj = lep.X()*lj_x+lep.Y()*lj_y+lep.Z()*lj_z;
                  // absolute value squared and normalized
                  float pLrel2 = lepXlj*lepXlj/lj2;
                  // lep2 = pTrel2 + pLrel2
                  pTrel2 = lep2-pLrel2;
                } else
                  std::cout << "JetCorrectionUncertainty not positive lepton-jet momentum: "<<lj2 << std::endl;
                return (pTrel2 > 0) ? std::sqrt(pTrel2) : 0.0;
              }


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

