/*! Definition of PatVertex class which contains analysis-level definiton of a vertex.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "DataFormats/VertexReco/interface/Vertex.h"

namespace pat {
class Vertex : public reco::Vertex {
public:
    typedef math::XYZVector Vector;

    Vertex() {}
    Vertex(const reco::Vertex& aVertex) : reco::Vertex(aVertex)
    {
        if(hasRefittedTracks()) {
            for(std::vector<reco::Track>::const_iterator iter = refittedTracks().begin();
                iter != refittedTracks().end(); ++iter) {
                trackMomentums.push_back(iter->momentum());
                trackWeights.push_back(trackWeight(originalTrack(*iter)));
            }
        }
        else {
            for(std::vector<reco::TrackBaseRef>::const_iterator iter = tracks_begin(); iter != tracks_end(); iter++) {
                trackMomentums.push_back((*iter)->momentum());
                trackWeights.push_back(trackWeight(*iter));
            }
        }
        _sumPt = __momentum(0.13957018, -1).pt();
        _sumPtSquared = __sumPtSquared(-1);
    }

    double sumPt() const { return _sumPt; }
    double sumPtSquared() const { return _sumPtSquared; }

    double __sumPtSquared(float minWeight = 0.5) const
    {
        double sum = 0.0;
        for(size_t n = 0; n < trackMomentums.size(); ++n) {
            if(trackWeights.at(n) >= minWeight)
                sum += std::sqrt( trackMomentums.at(n).Perp2() );
        }
        return sum;
    }

    math::XYZTLorentzVectorD __momentum(float mass = 0.13957018, float minWeight = 0.5) const
    {
        math::XYZTLorentzVectorD sum;
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double> > vec;

        for(size_t n = 0; n < trackMomentums.size(); ++n) {
            if(trackWeights.at(n) >= minWeight) {
                vec.SetPx(trackMomentums.at(n).x());
                vec.SetPy(trackMomentums.at(n).y());
                vec.SetPz(trackMomentums.at(n).z());
                vec.SetM(mass);
                sum += vec;
            }
        }
        return sum;
    }

    size_t numberOfTracks() const { return trackMomentums.size(); }
    size_t numberOfWeights() const { return trackWeights.size(); }

private:
    std::vector<Vector> trackMomentums;
    std::vector<float> trackWeights;
    double _sumPt;
    double _sumPtSquared;
};

typedef std::vector<Vertex> VertexCollection;

} // pat
