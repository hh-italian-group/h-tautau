#include "PhysicsTools/JetMCUtils/interface/JetMCTag.h"

#include "DataFormats/TauReco/interface/PFTauDecayMode.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TMath.h"
#include "TLorentzVector.h"


inline void findDaughters(const reco::GenParticle* mother, std::vector<const reco::GenParticle*>& daughters, int status)
{
  unsigned numDaughters = mother->numberOfDaughters();
  for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
    const reco::GenParticle* daughter = mother->daughterRef(iDaughter).get();

    if ( status == -1 || daughter->status() == status ) daughters.push_back(daughter);

    findDaughters(daughter, daughters, status);
  }
}

inline bool isNeutrino(const reco::GenParticle* daughter)
{
  return ( TMath::Abs(daughter->pdgId()) == 12 || TMath::Abs(daughter->pdgId()) == 14 || TMath::Abs(daughter->pdgId()) == 16 );
}

namespace utils_genMatch{
inline reco::Candidate::LorentzVector getVisMomentumNoLep(const std::vector<const reco::GenParticle*>& daughters, int status)
{
  reco::Candidate::LorentzVector p4Vis(0,0,0,0);

  for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
	daughter != daughters.end(); ++daughter ) {
    if (status != -1 && (*daughter)->status() != status) continue;
    if (isNeutrino(*daughter)) continue;
    if (TMath::Abs((*daughter)->pdgId()) == 11) continue;
    if (TMath::Abs((*daughter)->pdgId()) == 13) continue;
    //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
    //	  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
    p4Vis += (*daughter)->p4();
  }

  //std::cout << "--> vis. Momentum: Pt = " << p4Vis.pt() << ", eta = " << p4Vis.eta() << ", phi = " << p4Vis.phi() << std::endl;

  return p4Vis;
}

inline reco::Candidate::LorentzVector getVisMomentumNoLep(const reco::GenParticle* genLeg)
{
  std::vector<const reco::GenParticle*> stableDaughters;
  findDaughters(genLeg, stableDaughters, 1);

  reco::Candidate::LorentzVector p4Vis = getVisMomentumNoLep(stableDaughters,1);

  return p4Vis;
}
} //utils_genMatch

// inline const reco::GenParticle* findGenParticle(const reco::Candidate::LorentzVector& direction,
// 					 const reco::GenParticleCollection& genParticles, double dRmax, int status,
// 					 const std::vector<int>* pdgIds, bool pdgIdStrict)
// {
//   bool bestMatchMatchesPdgId = false;
//   const reco::GenParticle* bestMatch = 0;
//   
//   for ( reco::GenParticleCollection::const_iterator genParticle = genParticles.begin();
// 	genParticle != genParticles.end(); ++genParticle ) {
//     bool matchesPdgId = false;
//     if ( pdgIds ) {
//       for ( std::vector<int>::const_iterator pdgId = pdgIds->begin(); pdgId != pdgIds->end(); ++pdgId ) {
// 	if ( genParticle->pdgId() == (*pdgId) ) {
// 	  matchesPdgId = true;
// 	  break;
//         }
//       }
//     }
//     
//     // If we require strict PDG id checking, skip it if it doesn't match
//     if ( pdgIds && !matchesPdgId && pdgIdStrict ) continue;
//     
//     // Check if status matches - if not, skip it.
//     bool statusMatch = (status == -1 || genParticle->status() == status);
//     if ( !statusMatch ) continue;
// 
//     double dR = reco::deltaR(direction, genParticle->p4());
//     if ( dR > dRmax ) continue;
// 
//     // Check if higher than current best match
//     bool higherEnergyThanBestMatch = (bestMatch) ? 
//       (genParticle->energy() > bestMatch->energy()) : true;
// 
//     // Check if old bestMatch was not a prefered ID and the new one is.
//     if ( bestMatchMatchesPdgId ) {
//       // If the old one matches, only use the new one if it is a better
//       // energy match
//       if ( matchesPdgId && higherEnergyThanBestMatch ) bestMatch = &(*genParticle);
//     } else {
//       // The old one doesn't match the prefferred list if it is either
//       // a better energy match or better pdgId match
//       if ( higherEnergyThanBestMatch || matchesPdgId ) {
// 	bestMatch = &(*genParticle);
// 	if ( matchesPdgId ) bestMatchMatchesPdgId = true;
//       }
//     }
//   }
//   
//   return bestMatch;
// }

// inline reco::Candidate::Point getDecayVertex(const reco::GenParticle* mother)
// {
//   std::vector<const reco::GenParticle*> stableDaughters;
//   findDaughters(mother, stableDaughters, 1);
//   // If no daughters, return production vertex
//   if ( !stableDaughters.size() ) {
//     return mother->vertex();
//   } else {
//     // Try to sure we aren't going through any intermediate decays
//     for ( std::vector<const reco::GenParticle*>::const_iterator dau = stableDaughters.begin();
// 	  dau != stableDaughters.end(); ++dau ) {
//       if ( (*dau)->mother() == mother ) return (*dau)->vertex();
//     }
//     // If everything had an intermediate decay, just take the first one. :/
//     return stableDaughters[0]->vertex();
//   }
// }




// inline reco::Candidate::LorentzVector getVisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)
// {
//   reco::Candidate::LorentzVector p4Vis(0,0,0,0);
// 
//   for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
// 	daughter != daughters.end(); ++daughter ) {
//     if ( (status == -1 || (*daughter)->status() == status) && !isNeutrino(*daughter) ) {
//       //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
//       //	  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
//       p4Vis += (*daughter)->p4();
//     }
//   }
// 
//   //std::cout << "--> vis. Momentum: Pt = " << p4Vis.pt() << ", eta = " << p4Vis.eta() << ", phi = " << p4Vis.phi() << std::endl;
// 
//   return p4Vis;
// }
// 
// inline reco::Candidate::LorentzVector getVisMomentum(const reco::GenParticle* genLeg)
// {
//   std::vector<const reco::GenParticle*> stableDaughters;
//   findDaughters(genLeg, stableDaughters, 1);
// 
//   reco::Candidate::LorentzVector p4Vis = getVisMomentum(stableDaughters);
// 
//   return p4Vis;
// }
// 
// inline reco::Candidate::LorentzVector getInvisMomentum(const std::vector<const reco::GenParticle*>& daughters, int status)
// {
//   reco::Candidate::LorentzVector p4Invis(0,0,0,0);
// 
//   for ( std::vector<const reco::GenParticle*>::const_iterator daughter = daughters.begin();
// 	daughter != daughters.end(); ++daughter ) {
//     if ( (status == -1 || (*daughter)->status() == status) && isNeutrino(*daughter) ) {
//       //std::cout << "adding daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","
//       //	  << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
//       p4Invis += (*daughter)->p4();
//     }
//   }
// 
//   //std::cout << "--> invis. Momentum: Pt = " << p4Invis.pt() << ", phi = " << p4Invis.phi() << std::endl;
// 
//   return p4Invis;
// }
// 
// inline reco::Candidate::LorentzVector getInvisMomentum(const reco::GenParticle* genLeg)
// {
//   std::vector<const reco::GenParticle*> stableDaughters;
//   findDaughters(genLeg, stableDaughters, 1);
// 
//   reco::Candidate::LorentzVector p4Invis = getInvisMomentum(stableDaughters);
// 
//   return p4Invis;
// }
// 
// inline void compX1X2byCollinearApprox(double& x1, double& x2, double pxLeg1, double pyLeg1, double pxLeg2, double pyLeg2, double pxMEt, double pyMEt)
// {
//   double x1_numerator = pxLeg1*pyLeg2 - pxLeg2*pyLeg1;
//   double x1_denominator = pyLeg2*(pxLeg1 + pxMEt) - pxLeg2*(pyLeg1 + pyMEt);
//   x1 = ( x1_denominator != 0. ) ? x1_numerator/x1_denominator : -1.;
//   //std::cout << "x1 = " << x1 << std::endl;
// 
//   double x2_numerator = x1_numerator;
//   double x2_denominator = pxLeg1*(pyLeg2 + pyMEt) - pyLeg1*(pxLeg2 + pxMEt);
//   x2 = ( x2_denominator != 0. ) ? x2_numerator/x2_denominator : -1.;
//   //std::cout << "x2 = " << x2 << std::endl;
// }
// 
// inline double getPhysX(double x, bool& isWithinPhysRange)
// {
//   double physX = x;
// 
//   isWithinPhysRange = true;
// 
//   if ( x < 0. ) {
//     physX = 0.;
//     isWithinPhysRange = false;
//   }
// 
//   if ( x > 1. ) {
//     physX = 1.;
//     isWithinPhysRange = false;
//   }
// 
//   return physX;
// }
// 
// inline reco::Candidate::LorentzVector boostToRestFrame(const reco::Candidate::LorentzVector& p4daughter,
// 						const reco::Candidate::LorentzVector& p4mother)
// {
//   TLorentzVector p4mother_tlorentzvector(p4mother.px(), p4mother.py(), p4mother.pz(), p4mother.energy());
// 
//   TVector3 boostVector = -(p4mother_tlorentzvector.BoostVector());
// 
//   TLorentzVector p4daughter_tlorentzvector(p4daughter.px(), p4daughter.py(), p4daughter.pz(), p4daughter.energy());
//   p4daughter_tlorentzvector.Boost(boostVector);
// 
//   return reco::Candidate::LorentzVector(p4daughter_tlorentzvector.Px(), p4daughter_tlorentzvector.Py(), p4daughter_tlorentzvector.Pz(),
// 					p4daughter_tlorentzvector.Energy());
// }
// 
// //
// //-------------------------------------------------------------------------------
// //
// 
// inline void countDecayProducts(const reco::GenParticle* genParticle,
// 			int& numElectrons, int& numElecNeutrinos, int& numMuons, int& numMuNeutrinos, 
// 			int& numChargedHadrons, int& numPi0s, int& numOtherNeutralHadrons, int& numPhotons)
// {
//   int absPdgId = TMath::Abs(genParticle->pdgId());
//   int status   = genParticle->status();
//   int charge   = genParticle->charge();
// 
//   if      ( absPdgId == 111 ) ++numPi0s;
//   else if ( status   ==   1 ) {
//     if      ( absPdgId == 11 ) ++numElectrons;
//     else if ( absPdgId == 12 ) ++numElecNeutrinos;
//     else if ( absPdgId == 13 ) ++numMuons;
//     else if ( absPdgId == 14 ) ++numMuNeutrinos;
//     else if ( absPdgId == 15 ) { 
//       edm::LogError ("countDecayProducts")
//         << "Found tau lepton with status code 1 !!";
//       return; 
//     }
//     else if ( absPdgId == 16 ) return; // no need to count tau neutrinos
//     else if ( absPdgId == 22 ) ++numPhotons;
//     else if ( charge   !=  0 ) ++numChargedHadrons;
//     else                       ++numOtherNeutralHadrons;
//   } else {
//     unsigned numDaughters = genParticle->numberOfDaughters();
//     for ( unsigned iDaughter = 0; iDaughter < numDaughters; ++iDaughter ) {
//       const reco::GenParticle* daughter = genParticle->daughterRef(iDaughter).get();
// 
//       countDecayProducts(daughter, 
// 			 numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
// 			 numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons);
//     }
//   }
// }
// 
// inline std::string getGenTauDecayMode(const reco::GenParticle* genParticle) 
// {
// //--- determine generator level tau decay mode
// //
// //    NOTE: 
// //        (1) function implements logic defined in PhysicsTools/JetMCUtils/src/JetMCTag::genTauDecayMode
// //            for different type of argument 
// //        (2) this implementation should be more robust to handle cases of tau --> tau + gamma radiation
// //
//   
//   //std::cout << "<getGenTauDecayMode>:" << std::endl;
// 
//   int numElectrons           = 0;
//   int numElecNeutrinos       = 0;
//   int numMuons               = 0;
//   int numMuNeutrinos         = 0; 
//   int numChargedHadrons      = 0;
//   int numPi0s                = 0; 
//   int numOtherNeutralHadrons = 0;
//   int numPhotons             = 0;
// 
//   countDecayProducts(genParticle,
// 		     numElectrons, numElecNeutrinos, numMuons, numMuNeutrinos,
// 		     numChargedHadrons, numPi0s, numOtherNeutralHadrons, numPhotons);
// 
//   if      ( numElectrons == 1 && numElecNeutrinos == 1 ) return std::string("electron");
//   else if ( numMuons     == 1 && numMuNeutrinos   == 1 ) return std::string("muon");
//   
//   switch ( numChargedHadrons ) {
//   case 1 : 
//     if ( numOtherNeutralHadrons != 0 ) return std::string("oneProngOther");
//     switch ( numPi0s ) {
//     case 0:
//       return std::string("oneProng0Pi0");
//     case 1:
//       return std::string("oneProng1Pi0");
//     case 2:
//       return std::string("oneProng2Pi0");
//     default:
//       return std::string("oneProngOther");
//     }
//   case 3 : 
//     if ( numOtherNeutralHadrons != 0 ) return std::string("threeProngOther");
//     switch ( numPi0s ) {
//     case 0:
//       return std::string("threeProng0Pi0");
//     case 1:
//       return std::string("threeProng1Pi0");
//     default:
//       return std::string("threeProngOther");
//     }
//   default:
//     return std::string("rare");
//   }
// }
// 
// //
// //-------------------------------------------------------------------------------
// //
// 
// inline TVector2 getDiTauBisectorDirection(const reco::Candidate::LorentzVector& leg1P4, const reco::Candidate::LorentzVector& leg2P4)
// {
//   double leg1CosPhi = TMath::Cos(leg1P4.phi());
//   double leg1SinPhi = TMath::Sin(leg1P4.phi());
// 
//   double leg2CosPhi = TMath::Cos(leg2P4.phi());
//   double leg2SinPhi = TMath::Sin(leg2P4.phi());
// 
//   double diTauPx_unnormalized = leg1CosPhi + leg2CosPhi;
//   double diTauPy_unnormalized = leg1SinPhi + leg2SinPhi;
//   double diTauP = TMath::Sqrt(diTauPx_unnormalized*diTauPx_unnormalized + diTauPy_unnormalized*diTauPy_unnormalized);
//   double diTauPx_normalized = ( diTauP > 0. ) ? (diTauPx_unnormalized/diTauP) : diTauPx_unnormalized;
//   double diTauPy_normalized = ( diTauP > 0. ) ? (diTauPy_unnormalized/diTauP) : diTauPy_unnormalized;
// 
//   return TVector2(diTauPx_normalized, diTauPy_normalized);
// }
// 
// //
// //-------------------------------------------------------------------------------
// //
// 
// inline std::string getTauDecayModeName(int tauDecayMode)
// {
//   // "translate" from enum defined in DataFormats/TauReco/interface/PFTauDecayMode.h
//   // to string ( in format defined in PhysicsTools/JetMCUtils/src/JetMCTag.cc )
//   //
//   // NOTE: the "unphysical" 2-prong modes take into account
//   //       track reconstruction inefficiencies (migration from 3-prong decays),
//   //       fake tracks and tracks from the underlying event (migration from 1-prong decays)
//   //
//   if      ( tauDecayMode == reco::PFTauDecayMode::tauDecaysElectron           ) return std::string("electron");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecayMuon                ) return std::string("muon");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion0PiZero ) return std::string("oneProng0Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion1PiZero ) return std::string("oneProng1Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion2PiZero ) return std::string("oneProng2Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion3PiZero ) return std::string("oneProng3Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay1ChargedPion4PiZero ) return std::string("oneProng4Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion0PiZero ) return std::string("twoProng0Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion1PiZero ) return std::string("twoProng1Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion2PiZero ) return std::string("twoProng2Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion3PiZero ) return std::string("twoProng3Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay2ChargedPion4PiZero ) return std::string("twoProng4Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion0PiZero ) return std::string("threeProng0Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion1PiZero ) return std::string("threeProng1Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion2PiZero ) return std::string("threeProng2Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion3PiZero ) return std::string("threeProng3Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecay3ChargedPion4PiZero ) return std::string("threeProng4Pi0");
//   else if ( tauDecayMode == reco::PFTauDecayMode::tauDecayOther               ) return std::string("other");
//   else {
//     edm::LogError ("getTauDecayModeName")
//       << " Invalid tau decay mode = " << tauDecayMode << " !!";
//     return std::string("unknown");
//   }
// }
// 
// //
// //-------------------------------------------------------------------------------
// //
// 
// inline const reco::Candidate* getDistPion(const pat::Tau& recTauJet)
// {
//   if ( !recTauJet.isPFTau() ) {
//     //edm::LogWarning ("getDistPion")
//     //  << " Cannot identify 'distinguishable' pion for CaloTaus/TCTaus --> returning NULL pointer !!";
//     return 0;
//   }
// 
//   std::vector<reco::PFCandidatePtr> recTauJetChargedConstituents = recTauJet.signalPFChargedHadrCands();
// 
//   if ( recTauJetChargedConstituents.size() == 1 ) {
// 
// //--- tau- --> one-prong case (in particular rho- --> pi- pi0 or tau- --> a1- --> pi- pi0 pi0);
// //    the "distinguishable" pion is the leading charged hadron
//     return recTauJet.leadPFChargedHadrCand().get();
//   } else if ( recTauJetChargedConstituents.size() == 3 ) {
//     double recTauJetCharge = recTauJet.charge();
// 
//     for ( std::vector<reco::PFCandidatePtr>::const_iterator recTauJetChargedConstituent = recTauJetChargedConstituents.begin();
// 	  recTauJetChargedConstituent != recTauJetChargedConstituents.end(); ++recTauJetChargedConstituent ) {
// //--- tau- --> three-prong case (in particular a1- --> pi- pi+ pi-);
// //    the "distinguishable" pion is the pion of charge opposite to the tau-jet charge
//       if ( (*recTauJetChargedConstituent)->charge()*recTauJetCharge < 0 ) {
//         return recTauJetChargedConstituent->get();
//       }
//     }
//   } else {
//     //edm::LogWarning ("getDistPion")
//     //  << " Unsupported rec. tau decay mode = " << recTauJet.decayMode() << " --> returning NULL pointer !!";
//     return 0;
//   }
//   
//   //edm::LogWarning ("getDistPion")
//   //  << " Failed to identify 'distinguishable' rec. pion --> returning NULL pointer !!";
//   return 0;
// }
// 
// inline const reco::Candidate* getDistPion(const reco::GenJet& genTauJet)
// {
//   std::string genTauDecayMode = JetMCTagUtils::genTauDecayMode(genTauJet);
// 
//   std::vector<const reco::GenParticle*> genTauJetConstituents = genTauJet.getGenConstituents();
// 
//   if ( genTauDecayMode == "oneProng1Pi0" ||
//        genTauDecayMode == "oneProng2Pi0" ) {
//     for ( std::vector<const reco::GenParticle*>::const_iterator genTauJetConstituent = genTauJetConstituents.begin();
// 	  genTauJetConstituent != genTauJetConstituents.end(); ++genTauJetConstituent ) {
// 
// //--- tau- --> rho- --> pi- pi0 or tau- --> a1- --> pi- pi0 pi0 case;
// //    the "distinguishable" pion is the charged one
//       if ((*genTauJetConstituent)->charge()) return (*genTauJetConstituent);
//     }
//   } else if ( genTauDecayMode == "threeProng0Pi0" ) {
//     double genTauJetCharge = genTauJet.charge();
// 
//     for ( std::vector<const reco::GenParticle*>::const_iterator genTauJetConstituent = genTauJetConstituents.begin();
// 	  genTauJetConstituent != genTauJetConstituents.end(); ++genTauJetConstituent ) {
// 
// //--- tau- --> a1- --> pi- pi+ pi- case
// //    the "distinguishable" pion is the pion of charge opposite to the tau-jet charge
//       if ( (*genTauJetConstituent)->charge() != 0. &&
// 	   ((*genTauJetConstituent)->charge()*genTauJetCharge) < 0. ) {
// 	return *genTauJetConstituent;
//       }
//     }
//   } else {
//     //edm::LogWarning ("getDistPion")
//     //  << " Unsupported gen. tau decay mode = " << genTauDecayMode << " --> returning NULL pointer !!";
//     return 0;
//   }
//   //edm::LogWarning ("getDistPion")
//   //  << " Failed to identify 'distinguishable' gen. pion --> returning NULL pointer !!";
//   return 0;
// }
// 
// //
// //-------------------------------------------------------------------------------
// //
// 
// inline std::pair<double, double> compMEtProjU(const reco::Candidate::LorentzVector& zP4, double metPx, double metPy, int& errorFlag, bool subtract_qT)
// {
//   double qX = zP4.px();
//   double qY = zP4.py();
//   
//   return compMEtProjU(qX, qY, metPx, metPy, errorFlag, subtract_qT);
// }
// 
// inline std::pair<double, double> compMEtProjU(double qX, double qY, double metPx, double metPy, int& errorFlag, bool subtract_qT)
// {
//   double qT = TMath::Sqrt(qX*qX + qY*qY);
//   if ( qT == 0. ) {
//     edm::LogWarning ("compMEtProjU")
//       << " Failed to compute projection, because Z0 candidate has zero Pt --> returning dummy solution !!";
//     errorFlag = 1;
//     return std::pair<double, double>(0., 0.);
//   }
//   
//   double uX = -metPx;
//   double uY = -metPy;
//   if ( subtract_qT ) {
//     uX -= qX;
//     uY -= qY;
//   }
//   
//   double u1 = (uX*qX + uY*qY)/qT;
//   double u2 = (uX*qY - uY*qX)/qT;
//   
//   return std::pair<double, double>(u1, u2);
// }
// 
// //
// //-----------------------------------------------------------------------------------------------------------------------
// //
// 
// inline std::vector<double> compTrackPtSums(const reco::VertexCollection& vertices)
// {
//   size_t numVertices = vertices.size();
//   std::vector<double> trackPtSums(numVertices);
//   for ( size_t iVertex = 0; iVertex < numVertices; ++iVertex ) {
//     const reco::Vertex& vertex = vertices[iVertex];
// 
//     double trackPtSum = 0.;
//     for ( reco::Vertex::trackRef_iterator track = vertex.tracks_begin();
// 	  track != vertex.tracks_end(); ++track ) {
//       trackPtSum += (*track)->pt();
//     }
//     
//     trackPtSums[iVertex] = trackPtSum;
//   }
// 
//   return trackPtSums;
// }
// 
// inline size_t getNumVerticesPtGtThreshold(const std::vector<double>& trackPtSums, double ptThreshold)
// {
//   size_t numVertices = 0;
//   for ( std::vector<double>::const_iterator trackPtSum = trackPtSums.begin();
// 	trackPtSum != trackPtSums.end(); ++trackPtSum ) {
//     if ( (*trackPtSum) > ptThreshold ) ++numVertices;
//   }
//   
//   return numVertices;
// }
// 
// 
// inline reco::Candidate::LorentzVector getLeadChHadMomentum(const reco::GenParticle* genLeg)
// {
//   std::vector<const reco::GenParticle*> stableDaughters;
//   findDaughters(genLeg, stableDaughters, 1);
// 
//   reco::Candidate::LorentzVector p4LCH(0,0,0,0);
//   double lchPt = 0;
//   for ( std::vector<const reco::GenParticle*>::const_iterator daughter = stableDaughters.begin();
//         daughter != stableDaughters.end(); ++daughter ) {
//     if ( (*daughter)->status() == 1 && (*daughter)->charge() != 0 ) {
//       //std::cout << "daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","                                                                 
//       //          << " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
//       if((*daughter)->pt() > lchPt){
// 	lchPt = (*daughter)->pt();
// 	p4LCH = (*daughter)->p4();
//       }
//     }
//   }
// 
//   //std::cout << "--> lch Momentum: Pt = " << p4LCH.pt() << ", eta = " << p4LCH.eta() << ", phi = " << p4LCH.phi() << std::endl;                                                 
// 
//   return p4LCH;
// }