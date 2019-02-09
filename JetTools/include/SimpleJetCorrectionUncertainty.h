#pragma once

#include "JetCorrectorParameters.h"

namespace jec {

class SimpleJetCorrectionUncertainty {
public:
    SimpleJetCorrectionUncertainty();
    SimpleJetCorrectionUncertainty(const std::string& fDataFile);
    SimpleJetCorrectionUncertainty(const JetCorrectorParameters& fParameters);
    ~SimpleJetCorrectionUncertainty();
    SimpleJetCorrectionUncertainty(const SimpleJetCorrectionUncertainty&) = delete;
    SimpleJetCorrectionUncertainty& operator= (const SimpleJetCorrectionUncertainty&) = delete;

    const JetCorrectorParameters& parameters() const;
    float uncertainty(const std::vector<float>& fX, float fY, bool fDirection) const;

private:
    int findBin(const std::vector<float>& v, float x) const;
    float uncertaintyBin(unsigned fBin, float fY, bool fDirection) const;
    float linearInterpolation(float fZ, const float fX[2], const float fY[2]) const;

    JetCorrectorParameters* mParameters;
};

} // namespace jec
