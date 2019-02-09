#include "h-tautau/JetTools/include/SimpleJetCorrectionUncertainty.h"

#include <iostream>
#include "AnalysisTools/Core/include/exception.h"

namespace jec {

SimpleJetCorrectionUncertainty::SimpleJetCorrectionUncertainty()
{
    mParameters = new JetCorrectorParameters();
}

SimpleJetCorrectionUncertainty::SimpleJetCorrectionUncertainty(const std::string& fDataFile)
{
    mParameters = new JetCorrectorParameters(fDataFile);
}

SimpleJetCorrectionUncertainty::SimpleJetCorrectionUncertainty(const JetCorrectorParameters& fParameters)
{
    mParameters = new JetCorrectorParameters(fParameters);
}

SimpleJetCorrectionUncertainty::~SimpleJetCorrectionUncertainty()
{
    delete mParameters;
}

const JetCorrectorParameters& SimpleJetCorrectionUncertainty::parameters() const
{
    return *mParameters;
}

float SimpleJetCorrectionUncertainty::uncertainty(const std::vector<float>& fX, float fY, bool fDirection) const
{
    float result = 1.;
    int bin = mParameters->binIndex(fX);
    if (bin<0) {
        std::cout << "SimpleJetCorrectionUncertainty bin variables out of range" << std::endl;
        result = -999.0;
    } else
        result = uncertaintyBin(unsigned(bin),fY,fDirection);
    return result;
}

int SimpleJetCorrectionUncertainty::findBin(const std::vector<float>& v, float x) const
{
    size_t i;
    size_t n = v.size()-1;
    if (n<=0) return -1;
    if (x<v[0] || x>=v[n])
        return -1;
    for(i=0;i<n;i++) {
        if (x>=v[size_t(i)] && x<v[size_t(i+1)])
        return int(i);
    }
    return 0;
}

float SimpleJetCorrectionUncertainty::uncertaintyBin(unsigned fBin, float fY, bool fDirection) const
{
    if (fBin >= mParameters->size()) {
        std::cout << "SimpleJetCorrectionUncertainty:  wrong bin: "<<fBin<<": only "<<mParameters->size()
                  << " are available" << std::endl;
        return -999.0;
    }
    const std::vector<float>& p = mParameters->record(fBin).parameters();
    if ((p.size() % 3) != 0)
        throw analysis::exception("SimpleJetCorrectionUncertainty, wrong # of parameters: multiple of 3 expected,"
                                  " '%1%' got") % p.size();
    std::vector<float> yGrid,value;
    unsigned int N = unsigned(p.size()/3);
    float result = -1.0;
    for(unsigned i=0;i<N;i++) {
        unsigned ind = 3*i;
        yGrid.push_back(p[ind]);
        if (fDirection)// true = UP
            value.push_back(p[ind+1]);
        else // false = DOWN
            value.push_back(p[ind+2]);
    }
    if (fY <= yGrid[0])
        result = value[0];
    else if (fY >= yGrid[N-1])
        result = value[N-1];
    else {
        size_t bin = size_t(findBin(yGrid,fY));
        float vx[2],vy[2];
        for(size_t i=0;i<2;i++) {
            vx[i] = yGrid[bin+i];
            vy[i] = value[bin+i];
        }
        result = linearInterpolation(fY,vx,vy);
    }
    return result;
}

float SimpleJetCorrectionUncertainty::linearInterpolation(float fZ, const float fX[2], const float fY[2]) const
{
    // Linear interpolation through the points (x[i],y[i]). First find the line that
    // is defined by the points and then calculate the y(z).
    float r = 0;
    if (fX[0] == fX[1]) {
        if (fY[0] == fY[1])
            r = fY[0];
        else {
            std::cout << "SimpleJetCorrectionUncertainty: interpolation error" << std::endl;
            return -999.0;
        }
    } else {
        float a = (fY[1]-fY[0])/(fX[1]-fX[0]);
        float b = (fY[0]*fX[1]-fY[1]*fX[0])/(fX[1]-fX[0]);
        r = a*fZ+b;
    }
    return r;
}

} // namespace jec
