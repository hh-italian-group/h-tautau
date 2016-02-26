/*!
 * \file SkimFilter.cc
 * \brief Definition of SkimFilter class which applies skim for X->HH->bbTauTau PAT event selection.
 * \author Rosamaria Venditti (INFN Bari, Bari University)
 * \author Contributing author: Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Contributing author: Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-04-30 created
 *
 * Copyright 2014 Rosamaria Venditti,
 *                Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
//#include "HHbbTauTau/PatProduction/interface/PatVertex.h"

class SkimFilterMiniAOD : public edm::EDFilter {
public:
    explicit SkimFilterMiniAOD(const edm::ParameterSet& iConfig)
        : vertexTag(mayConsume<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("vertexSrc"))),
          muonTag(mayConsume<pat::MuonCollection >(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc"))),
          electronTag(mayConsume<pat::ElectronCollection>(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))),
          tauTag(mayConsume<pat::TauCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))) {}

private:
    virtual bool filter(edm::Event&, const edm::EventSetup&);

private:
//    const edm::EDGetTokenT<reco::VertexCollection> vertexTag;
    edm::EDGetToken vertexTag;
    edm::EDGetToken muonTag;
    edm::EDGetToken electronTag;
    edm::EDGetToken tauTag;
};

bool SkimFilterMiniAOD::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<reco::VertexCollection> hVertices;
    iEvent.getByToken(vertexTag, hVertices);
    const reco::VertexCollection& vertices = *hVertices;

    edm::Handle<pat::MuonCollection> hMuons;
    iEvent.getByToken(muonTag, hMuons);
    const pat::MuonCollection& muons = *hMuons.product();

    edm::Handle<pat::ElectronCollection> hElectrons;
    iEvent.getByToken(electronTag, hElectrons);
    const pat::ElectronCollection& electrons = *hElectrons.product();

    edm::Handle<pat::TauCollection> hTaus;
    iEvent.getByToken(tauTag, hTaus);
    const pat::TauCollection& taus = *hTaus.product();

    bool haveGoodVertex = false;
    unsigned nGoodTaus = 0;

  //  auto vtx = vertices.at(0);

    for (const reco::Vertex& vertex : vertices){
        if(vertex.ndof() > 4 && std::abs(vertex.z()) < 24 && vertex.position().rho() < 2){
            haveGoodVertex = true;
            break;
        }
    }
    if (!haveGoodVertex) {
        return false;}
  

    for(const pat::Tau& tau : taus) {
        if(tau.pt() > 18 && std::abs(tau.eta()) < 2.3  /*&& tau.tauID("decayModeFinding") > 0.5*/)
            ++nGoodTaus;
        if(nGoodTaus >= 2) return true;
    }
    if(!nGoodTaus) return false;
  

    for(const pat::Electron& electron : electrons) {
        if(electron.pt() > 20.0 && std::abs(electron.eta()) < 2.0) {
            std::cout<<" Ele Missing Hits =  "<<
                       electron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS)<<std::endl;
            return true;}
    }

    for(const pat::Muon& muon : muons) {
//        bool muonIP = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.045 &&
//                      fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.2;
        if(muon.pt() > 16.0 && std::abs(muon.eta())<2.4 && muon.isLooseMuon()) {
//            float iso = (muon.pfIsolationR03().sumChargedHadronPt + std::max(
//                       muon.pfIsolationR03().sumNeutralHadronEt +
//                       muon.pfIsolationR03().sumPhotonEt -
//                       0.5 * muon.pfIsolationR03().sumPUPt, 0.0)) / muon.pt();

//            std::cout<<" Muon Isolament =  "<<iso<<std::endl;

            return true;
        }
    }

    return false;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SkimFilterMiniAOD);
