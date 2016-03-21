/*! Definition of candidate class to store reconstructed object candidate.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <vector>

#include <TLorentzVector.h>

#include "Particles.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"


namespace ntuple {
namespace tau_id {
enum hadronicDecayMode {
  kNull = -1,
  kOneProng0PiZero,
  kOneProng1PiZero,
  kOneProng2PiZero,
  kOneProng3PiZero,
  kOneProngNPiZero,
  kTwoProng0PiZero,
  kTwoProng1PiZero,
  kTwoProng2PiZero,
  kTwoProng3PiZero,
  kTwoProngNPiZero,
  kThreeProng0PiZero,
  kThreeProng1PiZero,
  kThreeProng2PiZero,
  kThreeProng3PiZero,
  kThreeProngNPiZero,
  kRareDecayMode
};

template<typename Value>
inline hadronicDecayMode ConvertToHadronicDecayMode(const Value& value)
{
    if(value < kNull || value > kRareDecayMode)
        throw std::runtime_error("value is not a hadronicDecayMode");
    return static_cast<ntuple::tau_id::hadronicDecayMode>(value);
}

}
}

namespace analysis {

//template<typename NtupleObject>
class CandidateV2;
typedef std::shared_ptr<const CandidateV2> CandidateV2Ptr;
typedef std::vector<CandidateV2Ptr> CandidateV2PtrVector;

//template<typename NtupleObject>
class CandidateV2 {
public:
    enum class Type { Electron, Muon, Tau, Jet, Z, Higgs };
    static int UnknownCharge() { return std::numeric_limits<int>::max(); }

private:
    template<typename PATObject>
    static Type TypeFromNtupleObject(const PATObject& ntupleObject)
    {
        static const std::map<size_t, Type> ntupleTypes = {
            { typeid(pat::Electron).hash_code(), Type::Electron }, { typeid(pat::Muon).hash_code(), Type::Muon },
            { typeid(pat::Tau).hash_code(), Type::Tau }, { typeid(pat::Jet).hash_code(), Type::Jet }
        };

        const std::type_info& objectType = typeid(ntupleObject);
        if(!ntupleTypes.count(objectType.hash_code()))
            throw exception("Unknown n-tuple type '%1%'.") % objectType.name();
        return ntupleTypes.at(objectType.hash_code());
    }

    Type type;
    int charge;
    TLorentzVector momentum;
    bool has_vertexPosition;
    TVector3 vertexPosition;

    CandidateV2PtrVector daughters;
    CandidateV2PtrVector finalStateDaughters;

    const reco::LeafCandidate* ntupleObject;

public:
    template<typename PATObject>
    explicit CandidateV2(const PATObject& _ntupleObject) :
        type(TypeFromNtupleObject<PATObject>(_ntupleObject)), charge(_ntupleObject.charge()),
        momentum(MakeLorentzVectorPtEtaPhiM(_ntupleObject.pt(), _ntupleObject.eta(), _ntupleObject.phi(), _ntupleObject.mass())),
        has_vertexPosition(true), vertexPosition(_ntupleObject.vx(), _ntupleObject.vy(), _ntupleObject.vz()),
        ntupleObject(&_ntupleObject)
    {}

    explicit CandidateV2(const pat::Jet& _ntupleObject) :
        type(TypeFromNtupleObject<pat::Jet>(_ntupleObject)), charge(UnknownCharge()),
        momentum(MakeLorentzVectorPtEtaPhiM(_ntupleObject.pt(), _ntupleObject.eta(), _ntupleObject.phi(), _ntupleObject.mass())),
        has_vertexPosition(false), ntupleObject(&_ntupleObject) {}

    CandidateV2(Type _type, const CandidateV2Ptr& daughter1, const CandidateV2Ptr& daughter2)
        : type(_type), has_vertexPosition(false), ntupleObject(nullptr)
    {
        if(!daughter1 || !daughter2)
            throw exception("Candidate daughters can't be nullptr.");
        if(daughter1 == daughter2)
            throw exception("Can't create composit candidate use identical daughters.");

        daughters.push_back(daughter1);
        daughters.push_back(daughter2);
        momentum = daughter1->momentum + daughter2->momentum;
        if (daughter1->charge == UnknownCharge() || daughter2->charge == UnknownCharge())
            charge = UnknownCharge();
        else charge = daughter1->charge + daughter2->charge;
        for(const auto& daughter : daughters) {
            if(!daughter->finalStateDaughters.size())
                finalStateDaughters.push_back(daughter);
            for(const auto& finalStateDaughter : daughter->finalStateDaughters)
                finalStateDaughters.push_back(finalStateDaughter);
        }
    }

    Type GetType() const { return type; }
    int GetCharge() const { return charge; }
    const TLorentzVector& GetMomentum() const { return momentum; }
    bool HasVertexPosition() const { return has_vertexPosition; }
    const TVector3& GetVertexPosition() const
    {
        if(!HasVertexPosition())
            throw exception("Vertex position is not defined.");
        return vertexPosition;
    }

    template<typename PATObject>
    const PATObject& GetNtupleObject() const
    {
        if(!ntupleObject)
            throw exception("Candidate is not associated with antuple object.");
        const PATObject* casted = dynamic_cast<const PATObject*>(ntupleObject);
        if(!casted)
            throw exception("Bad ntuple object type '%1%'. Expected '%2%'.") % typeid(PATObject).name()
                                                                             % typeid(ntupleObject).name();
        return *casted;
    }

    const CandidateV2PtrVector& GetDaughters() const { return daughters; }
    const CandidateV2PtrVector& GetFinalStateDaughters() const { return finalStateDaughters; }

    CandidateV2Ptr GetDaughter(Type daughterType) const
    {
        const auto selected_daughters = GetDaughters(daughterType);
        if(!selected_daughters.size())
            throw exception("Daughter with specified type not found.");
        if(selected_daughters.size() > 1)
            throw exception("More than one daughter with specified type are found.");
        return selected_daughters.front();
    }

    CandidateV2PtrVector GetDaughters(Type daughterType) const
    {
        CandidateV2PtrVector result;
        for(const auto& daughter : daughters) {
            if(daughter->type == daughterType)
                result.push_back(daughter);
        }
        return result;
    }

    const CandidateV2Ptr& GetLeadingDaughter(Type expectedDaughterType) const
    {
        const auto& leadingDaughter = daughters.at(GetPtOrderedDaughterIndexes().first);
        if(leadingDaughter->type != expectedDaughterType)
            throw exception("Unexpected leading daughter type.");
        return leadingDaughter;
    }

    const CandidateV2Ptr& GetSubleadingDaughter(Type expectedDaughterType) const
    {
        const auto& subleadingDaughter = daughters.at(GetPtOrderedDaughterIndexes().second);
        if(subleadingDaughter->type != expectedDaughterType)
            throw exception("Unexpected subleading daughter type.");
        return subleadingDaughter;
    }

    bool operator != (const CandidateV2& other) const
    {
        if (!daughters.size() && !other.daughters.size()) return this->type != other.type;
        if (daughters.size() != other.daughters.size()) return true;
        for (size_t n = 0; n < daughters.size(); ++n){
            CandidateV2Ptr daughter = daughters.at(n);
            CandidateV2Ptr other_daughter = other.daughters.at(n);
            if (*daughter != *other_daughter) return true;
        }
        return false;
    }

    bool operator == (const CandidateV2& other) const
    {
        return !(*this != other);
    }

private:
    std::pair<size_t, size_t> GetPtOrderedDaughterIndexes() const
    {
        if(!daughters.size())
            throw exception("Candidate has no daughters.");
        return daughters.at(0)->momentum.Pt() >= daughters.at(1)->momentum.Pt()
                ? std::pair<size_t, size_t>(0, 1) : std::pair<size_t, size_t>(1, 0);
    }
};

class VertexV2;
typedef std::shared_ptr<const VertexV2> VertexV2Ptr;
typedef std::vector<VertexV2Ptr> VertexV2PtrVector;

class VertexV2 {
public:
    VertexV2(const reco::Vertex& _ntupleObject)
        : ndf(_ntupleObject.ndof()),
          position(_ntupleObject.x(), _ntupleObject.y(), _ntupleObject.z()), ntupleObject(&_ntupleObject) {}

    bool operator< (const VertexV2& other) const { return ntupleObject < other.ntupleObject; }

    //double GetSumPtSquared() const { return sumPtSquared; }
    unsigned GetNdf() const { return ndf; }
    const TVector3& GetPosition() const { return position; }
    const reco::Vertex& GetNtupleObject() const { return *ntupleObject; }

private:
   // double sumPtSquared;
    unsigned ndf;
    TVector3 position;
    const reco::Vertex* ntupleObject;
};

class MissingET;
typedef std::shared_ptr<const MissingET> MissingETPtr;

class MissingET {
public:
    template<typename METtype>
    MissingET(const METtype& _ntupleObject)
        :pt(_ntupleObject.pt()), phi(_ntupleObject.phi()), isPF(_ntupleObject.isPFMET()),
          ntupleObject(&_ntupleObject) {
        CovVector = METCovMatrixToVector(_ntupleObject.getSignificanceMatrix());
    }

    template<typename METtype>
    MissingET(const METtype& _ntupleObject,const reco::METCovMatrix& m)
        :pt(_ntupleObject.pt()), px(_ntupleObject.px()), py(_ntupleObject.py()),
         phi(_ntupleObject.phi()), isPF(_ntupleObject.isPFMET()), ntupleObject(&_ntupleObject) {
        CovVector = METCovMatrixToVector(m);
    }

    //const reco::RecoCandidate& GetNtupleObject() const { return *ntupleObject; }
    const double& Pt() const  { return pt; }
    const double& Px() const  { return px; }
    const double& Py() const  { return py; }
    const double& Phi() const { return phi; }
    const bool&   isPFMET() const { return isPF; }
    const std::vector<Float_t> GetCovVector() const {return CovVector;}

    template<typename METtype>
    const METtype& GetNtupleObject() const
    {
        if(!ntupleObject)
            throw exception("Candidate is not associated with antuple object.");
        const METtype* casted = dynamic_cast<const METtype*>(ntupleObject);
        if(!casted)
            throw exception("Bad ntuple object type '%1%'. Expected '%2%'.") % typeid(METtype).name()
                                                                             % typeid(ntupleObject).name();

        return *casted;
    }

private:
    double pt,px,py;
    double phi;
    bool isPF;
    const reco::RecoCandidate* ntupleObject;
    //Matrix
    std::vector<Float_t> CovVector;
private:
    std::vector<Float_t> METCovMatrixToVector (const reco::METCovMatrix& m){
        std::vector<Float_t> v(4);
        v[0] = m[0][0]; // xx
        //std::cout<<"\t xx = "<< m[0][0] << "\t xy = "<<m[0][1] << "\n\t yx = "<< m[1][0] << "\t yy = "<< m[1][1] <<std::endl;
        v[1] = m[0][1]; // xy
        v[2] = m[1][0]; // yx
        v[3] = m[1][1]; // yy
        return v;
    }
};

namespace detail {
const std::map<CandidateV2::Type, std::string> CandidateV2TypeNameMap = {
    { CandidateV2::Type::Electron, "electron" }, { CandidateV2::Type::Muon, "muon" }, { CandidateV2::Type::Tau, "tau" },
    { CandidateV2::Type::Jet, "jet" }, { CandidateV2::Type::Z, "Z" }, { CandidateV2::Type::Higgs, "Higgs" }
};
} // namespace detail

// enum class Type { Electron, Muon, Tau, Jet, Z, Higgs };

std::ostream& operator<< (std::ostream& s, const CandidateV2::Type& t)
{
    s << detail::CandidateV2TypeNameMap.at(t);
    return s;
}


} // analysis

