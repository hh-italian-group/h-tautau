/*! Definition of Candidate class to store reconstructed object candidate.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include <vector>

#include <TLorentzVector.h>

#include "EventDescriptor.h"
#include "Particles.h"
#include "AnalysisTools/Core/include/exception.h"
#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {

class Candidate;
typedef std::shared_ptr<const Candidate> CandidatePtr;
typedef std::vector<CandidatePtr> CandidatePtrVector;

class Candidate {
public:
    enum class Type { Electron, Muon, Tau, Jet, Z, Higgs };
    static int UnknownCharge() { return std::numeric_limits<int>::max(); }

private:
    static Type TypeFromNtupleObject(const root_ext::detail::BaseDataClass& ntupleObject)
    {
        static const std::map<size_t, Type> ntupleTypes = {
            { typeid(ntuple::Electron).hash_code(), Type::Electron }, { typeid(ntuple::Muon).hash_code(), Type::Muon },
            { typeid(ntuple::Tau).hash_code(), Type::Tau }, { typeid(ntuple::Jet).hash_code(), Type::Jet }
        };

        const std::type_info& objectType = typeid(ntupleObject);
        if(!ntupleTypes.count(objectType.hash_code()))
            throw exception("Unknown n-tuple type '") << objectType.name() << "'.";
        return ntupleTypes.at(objectType.hash_code());
    }

    Type type;
    int charge;
    TLorentzVector momentum;
    bool has_vertexPosition;
    TVector3 vertexPosition;

    CandidatePtrVector daughters;
    CandidatePtrVector finalStateDaughters;
    const root_ext::detail::BaseDataClass* ntupleObject;

public:
    template<typename NtupleObject>
    explicit Candidate(const NtupleObject& _ntupleObject) :
        type(TypeFromNtupleObject(_ntupleObject)), charge(_ntupleObject.charge),
        momentum(MakeLorentzVectorPtEtaPhiM(_ntupleObject.pt, _ntupleObject.eta, _ntupleObject.phi, _ntupleObject.mass)),
        has_vertexPosition(true), vertexPosition(_ntupleObject.vx, _ntupleObject.vy, _ntupleObject.vz),
        ntupleObject(&_ntupleObject) {}

    explicit Candidate(const ntuple::Jet& _ntupleObject) :
        type(TypeFromNtupleObject(_ntupleObject)), charge(UnknownCharge()),
        momentum(MakeLorentzVectorPtEtaPhiM(_ntupleObject.pt, _ntupleObject.eta, _ntupleObject.phi, _ntupleObject.mass)),
        has_vertexPosition(false), ntupleObject(&_ntupleObject) {}

    Candidate(Type _type, const CandidatePtr& daughter1, const CandidatePtr& daughter2)
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

    template<typename NtupleObject>
    const NtupleObject& GetNtupleObject() const
    {
        if(!ntupleObject)
            throw exception("Candidate is not associated with antuple object.");
        const NtupleObject* casted = dynamic_cast<const NtupleObject*>(ntupleObject);
        if(!casted)
            throw exception("Bad ntuple object type '") << typeid(NtupleObject).name() << "'. Expected '"
                                                        << typeid(ntupleObject).name() << "'.";
        return *casted;
    }

    const CandidatePtrVector& GetDaughters() const { return daughters; }
    const CandidatePtrVector& GetFinalStateDaughters() const { return finalStateDaughters; }

    CandidatePtr GetDaughter(Type daughterType) const
    {
        const auto selected_daughters = GetDaughters(daughterType);
        if(!selected_daughters.size())
            throw exception("Daughter with specified type not found.");
        if(selected_daughters.size() > 1)
            throw exception("More than one daughter with specified type are found.");
        return selected_daughters.front();
    }

    CandidatePtrVector GetDaughters(Type daughterType) const
    {
        CandidatePtrVector result;
        for(const auto& daughter : daughters) {
            if(daughter->type == daughterType)
                result.push_back(daughter);
        }
        return result;
    }

    const CandidatePtr& GetLeadingDaughter(Type expectedDaughterType) const
    {
        const auto& leadingDaughter = daughters.at(GetPtOrderedDaughterIndexes().first);
        if(leadingDaughter->type != expectedDaughterType)
            throw exception("Unexpected leading daughter type.");
        return leadingDaughter;
    }

    const CandidatePtr& GetSubleadingDaughter(Type expectedDaughterType) const
    {
        const auto& subleadingDaughter = daughters.at(GetPtOrderedDaughterIndexes().second);
        if(subleadingDaughter->type != expectedDaughterType)
            throw exception("Unexpected subleading daughter type.");
        return subleadingDaughter;
    }

    bool operator != (const Candidate& other) const
    {
        if (!daughters.size() && !other.daughters.size()) return ntupleObject != other.ntupleObject;
        if (daughters.size() != other.daughters.size()) return true;
        for (size_t n = 0; n < daughters.size(); ++n){
            CandidatePtr daughter = daughters.at(n);
            CandidatePtr other_daughter = other.daughters.at(n);
            if (*daughter != *other_daughter) return true;
        }
        return false;
    }

    bool operator == (const Candidate& other) const
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

class Vertex;
typedef std::shared_ptr<const Vertex> VertexPtr;
typedef std::vector<VertexPtr> VertexPtrVector;

class Vertex {
public:
    Vertex(const ntuple::Vertex& _ntupleObject)
        : sumPtSquared(_ntupleObject.sumPtSquared), ndf(_ntupleObject.ndf),
          position(_ntupleObject.x, _ntupleObject.y, _ntupleObject.z), ntupleObject(&_ntupleObject) {}

    bool operator< (const Vertex& other) const { return ntupleObject < other.ntupleObject; }

    double GetSumPtSquared() const { return sumPtSquared; }
    unsigned GetNdf() const { return ndf; }
    const TVector3& GetPosition() const { return position; }
    const ntuple::Vertex& GetNtupleObject() const { return *ntupleObject; }

private:
    double sumPtSquared;
    unsigned ndf;
    TVector3 position;
    const ntuple::Vertex* ntupleObject;
};

namespace detail {
const std::map<Candidate::Type, std::string> CandidateTypeNameMap = {
    { Candidate::Type::Electron, "electron" }, { Candidate::Type::Muon, "muon" }, { Candidate::Type::Tau, "tau" },
    { Candidate::Type::Jet, "jet" }, { Candidate::Type::Z, "Z" }, { Candidate::Type::Higgs, "Higgs" }
};
} // namespace detail

// enum class Type { Electron, Muon, Tau, Jet, Z, Higgs };

std::ostream& operator<< (std::ostream& s, const Candidate::Type& t)
{
    s << detail::CandidateTypeNameMap.at(t);
    return s;
}


} // analysis
