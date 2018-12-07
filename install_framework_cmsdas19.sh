#!/usr/bin/env bash
# Install hh-italian-group framework.
# This file is part of https://github.com/hh-italian-group/h-tautau.

MODE="cmsdas19"
N_JOBS=4
RELEASE="CMSSW_9_4_12"
export SCRAM_ARCH="slc6_amd64_gcc630"

scramv1 project CMSSW $RELEASE
RESULT=$?
if [ $RESULT -ne 0 ] ; then
    echo "ERROR: unable to create working area for CMSSW release '$RELEASE'."
    exit 2
fi

cd $RELEASE/src
eval `scramv1 runtime -sh`
RESULT=$?
if [ $RESULT -ne 0 ] ; then
    echo "ERROR: unable to setup the environment for CMSSW release '$RELEASE'."
    exit 2
fi

# SVfit packages
git clone https://github.com/hh-italian-group/SVfit_standalone.git TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
git checkout hh_italian
cd ../..

# HHKinFit2 packages
git clone https://github.com/hh-italian-group/HHKinFit2.git HHKinFit2/HHKinFit2

# LeptonEfficiencies packages
git clone https://github.com/hh-italian-group/LeptonEff-interface.git HTT-utilities
git clone https://github.com/hh-italian-group/LeptonEfficiencies.git HTT-utilities/LepEffInterface/data

# Recoil Corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections

# Install analysis packages
declare -A ANA_PACKAGES
ANA_PACKAGES=( ["AnalysisTools"]="cmsdas19:cmsdas_2019" \
               ["h-tautau"]="cmsdas19:cmsdas_2019" \
               ["hh-bbtautau"]="cmsdas19:cmsdas_2019" )

for pkg in "${!ANA_PACKAGES[@]}" ; do
    pkg_descs="${ANA_PACKAGES[$pkg]}"
    branch="master"
    for desc in $pkg_descs ; do
        if [ "${desc%%:*}" = "$MODE" ] ; then
            branch=${desc##*:}
            break
        fi
    done

    git clone https://github.com/hh-italian-group/${pkg}.git
    cd "$pkg"
    if [ "$branch" != "master" ] ; then
        git checkout -b $branch origin/$branch
    fi
    cd ..
done

BUILD_PATH="../build"

# Prepare analysis working area
./AnalysisTools/Run/install.sh "$BUILD_PATH" "${!ANA_PACKAGES[@]}"

scram b -j$N_JOBS
