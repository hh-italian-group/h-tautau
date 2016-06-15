/*! Definition of candidate class to store reconstructed object candidate.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {

class AnalysisObject {
public:
    static constexpr int UnknownCharge = std::numeric_limits<int>::max();

    AnalysisObject() {}
    template<typename FourVector>
    explicit AnalysisObject(const FourVector& _momentum, int _charge = UnknownCharge)
        : momentum(_momentum), charge(_charge) {}
    virtual ~AnalysisObject() {}

    const LorentzVector& GetMomentum() const { return momentum; }
    template<typename FourVector>
    void SetMomentum(const FourVector& _momentum) { momentum = LorentzVector(_momentum); }

    int GetCharge() const { return charge; }
    void SetCharge(int _charge) { charge = _charge; }
    bool HasCharge() const { return charge != UnknownCharge; }

    bool operator<(const AnalysisObject& other) const { return momentum.pt() > other.momentum.pt(); }

private:
    LorentzVector momentum;
    int charge;
};

template<typename _PATObject, typename _PATObjectConstPtr = const _PATObject*>
class Candidate : public AnalysisObject {
public:
    using PATObject = _PATObject;
    using PATObjectConstPtr = _PATObjectConstPtr;

    explicit Candidate(const PATObject& _patObject) : AnalysisObject(_patObject.p4()), patObject(&_patObject) {}
    explicit Candidate(const PATObjectConstPtr& _patObjectPtr)
        : AnalysisObject(_patObjectPtr->p4()), patObject(&(*_patObjectPtr)), patObjectPtr(_patObjectPtr) {}

    const PATObject& operator*() const { return *patObject; }
    const PATObject* operator->() const { return patObject; }
    const PATObjectConstPtr& getPtr() const { return patObjectPtr; }
    bool operator==(const Candidate<PATObject, PATObjectConstPtr>& other) const { return patObject == other.patObject; }
    bool operator!=(const Candidate<PATObject, PATObjectConstPtr>& other) const { return !(*this == other); }

private:
    const PATObject* patObject;
    PATObjectConstPtr patObjectPtr;
};

namespace detail {
template<typename PATObject>
bool CompareIsolations(double iso_1, double iso_2);
}

template<typename PATObject, typename PATObjectConstPtr = const PATObject*>
class LeptonCandidate : public Candidate<PATObject, PATObjectConstPtr> {
public:
    LeptonCandidate(const PATObject& _patObject, double _isolation)
        : Candidate<PATObject, PATObjectConstPtr>(_patObject), isolation(_isolation)
    {
        this->SetCharge(_patObject.charge());
    }

    LeptonCandidate(const PATObjectConstPtr& _patObjectPtr, double _isolation)
        : Candidate<PATObject, PATObjectConstPtr>(_patObjectPtr), isolation(_isolation)
    {
        this->SetCharge(_patObjectPtr->charge());
    }

    double GetIsolation() const { return isolation; }

    bool IsMoreIsolated(const LeptonCandidate<PATObject, PATObjectConstPtr>& other) const
    {
        return detail::CompareIsolations<PATObject>(isolation, other.isolation);
    }

private:
    double isolation;
};

template<typename _FirstDaughter, typename _SecondDaughter>
class CompositCandidate : public AnalysisObject {
public:
    using FirstDaughter = _FirstDaughter;
    using SecondDaughter = _SecondDaughter;
    CompositCandidate(const FirstDaughter& _firstDaughter, const SecondDaughter& _secondDaughter)
        : AnalysisObject(_firstDaughter.GetMomentum() + _secondDaughter.GetMomentum()),
          firstDaughter(_firstDaughter), secondDaughter(_secondDaughter)
    {
        if(firstDaughter.HasCharge() && secondDaughter.HasCharge())
            this->SetCharge(firstDaughter.GetCharge() + secondDaughter.GetCharge());
        daughterMomentums.push_back(firstDaughter.GetMomentum());
        daughterMomentums.push_back(secondDaughter.GetMomentum());
    }

    const FirstDaughter& GetFirstDaughter() const { return firstDaughter; }
    const SecondDaughter& GetSecondDaughter() const { return secondDaughter; }
    const std::vector<LorentzVector>& GetDaughterMomentums() const { return daughterMomentums; }

private:
    FirstDaughter firstDaughter;
    SecondDaughter secondDaughter;
    std::vector<LorentzVector> daughterMomentums;
};

template<typename MetObject>
class MissingET : public AnalysisObject {
public:
    using CovMatrix = SquareMatrix<2>;
    MissingET(const MetObject& _metObject, const CovMatrix& _covMatrix)
        : AnalysisObject(LorentzVectorM(_metObject.pt(), 0, _metObject.phi(), 0)),
          metObject(&_metObject), covMatrix(_covMatrix)
    {
    }

    const MetObject* operator*() const { return metObject; }
    const MetObject* operator->() const { return metObject; }
    const CovMatrix& GetCovMatrix() const { return covMatrix; }

private:
    const MetObject* metObject;
    CovMatrix covMatrix;
};

} // analysis
