/*! Definition of candidate class to store reconstructed object candidate.
This file is part of https://github.com/hh-italian-group/h-tautau. */

#pragma once

#include "AnalysisTools/Core/include/AnalysisMath.h"

namespace analysis {

class AnalysisObject {
public:
    static constexpr int UnknownCharge = std::numeric_limits<int>::max();
    static constexpr double UnknownIsolation = std::numeric_limits<double>::max();

    AnalysisObject() : charge(UnknownCharge), isolation(UnknownIsolation) {}
    template<typename FourVector>
    explicit AnalysisObject(const FourVector& _momentum, int _charge = UnknownCharge,
                            double _isolation = UnknownIsolation) :
        momentum(_momentum), charge(_charge), isolation(_isolation) {}
    AnalysisObject(const AnalysisObject& other) = default;
    virtual ~AnalysisObject() {}
    AnalysisObject& operator=(const AnalysisObject& other) = default;

    const LorentzVector& GetMomentum() const { return momentum; }
    template<typename FourVector>
    void SetMomentum(const FourVector& _momentum) { momentum = LorentzVector(_momentum); }

    int GetCharge() const { return charge; }
    void SetCharge(int _charge) { charge = _charge; }
    bool HasCharge() const { return charge != UnknownCharge; }

    double GetIsolation() const { return isolation; }
    void SetIsolation(double _isolation) { isolation = _isolation; }
    bool HasIsolation() const { return isolation != UnknownIsolation; }

    bool operator<(const AnalysisObject& other) const { return momentum.pt() > other.momentum.pt(); }
private:
    LorentzVector momentum;
    int charge;
    double isolation;
};

template<typename _PATObject, typename _PATObjectConstPtr = const _PATObject*>
class Candidate : public AnalysisObject {
public:
    using PATObject = _PATObject;
    using PATObjectConstPtr = _PATObjectConstPtr;

    explicit Candidate(const PATObject& _patObject, int _charge = UnknownCharge, double _isolation = UnknownIsolation) :
        AnalysisObject(_patObject.p4(), _charge, _isolation), patObject(&_patObject)
    {
    }

    explicit Candidate(const PATObjectConstPtr& _patObjectPtr, int _charge = UnknownCharge,
                       double _isolation = UnknownIsolation) :
        AnalysisObject(_patObjectPtr->p4(), _charge, _isolation), patObject(&(*_patObjectPtr)),
        patObjectPtr(_patObjectPtr)
    {
    }

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
        : Candidate<PATObject, PATObjectConstPtr>(_patObject, _patObject.charge(), _isolation)
    {
    }

    LeptonCandidate(const PATObjectConstPtr& _patObjectPtr, double _isolation)
        : Candidate<PATObject, PATObjectConstPtr>(_patObjectPtr, _patObjectPtr->charge(), _isolation)
    {
    }

    bool IsMoreIsolated(const LeptonCandidate<PATObject, PATObjectConstPtr>& other) const
    {
        return detail::CompareIsolations<PATObject>(this->GetIsolation(), other.GetIsolation());
    }
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

    template<typename FourVector>
    CompositCandidate(const FirstDaughter& _firstDaughter, const SecondDaughter& _secondDaughter,
                      const FourVector& _momentum)
        : AnalysisObject(_momentum), firstDaughter(_firstDaughter), secondDaughter(_secondDaughter)
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
