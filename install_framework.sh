#!/bin/bash
# Install hh-italian-group framework.
# This file is part of https://github.com/hh-italian-group/hh-bbtautau.

INSTALL_MODES=(prod prod17 ana limits)
DEFAULT_N_JOBS=4
DEFAULT_RELEASE_PROD="CMSSW_8_0_28"
DEFAULT_RELEASE_PROD17="CMSSW_9_4_4"
DEFAULT_RELEASE_LIMITS="CMSSW_7_4_7"
DEFAULT_RELEASE_ANA="CMSSW_9_0_0"

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

if [ $# -lt 1 -o $# -gt 3 ] ; then
    echo "Usage: mode [n_jobs] [cmssw_release]"
    printf "\n\tmode\t\t\tinstallation mode. Supported modes: "
    join_by ", " "${INSTALL_MODES[@]}"
    printf ".\n\tn_jobs\t\t\tthe number of jobs to run simultaneous during the compilation. Default: $DEFAULT_N_JOBS.\n"
    printf "\tcmssw_release\t\tCMSSW release."
    printf " Default: $DEFAULT_RELEASE_PROD\tfor tuple production,\n"
    printf "\t\t\t\t\t\t\t$DEFAULT_RELEASE_LIMITS\tfor limits computation,\n"
    printf "\t\t\t\t\t\t\t$DEFAULT_RELEASE_ANA\tfor analysis.\n"
    exit 1
fi

MODE=$1
if [[ ! ${INSTALL_MODES[@]} =~ $MODE ]] ; then
    echo "ERROR: unsupported installation mode '$MODE'."
    echo "Supported installation modes: ${INSTALL_MODES[@]}"
    exit 1
fi

N_JOBS=$2
if [ "x$N_JOBS" = "x" ] ; then N_JOBS=$DEFAULT_N_JOBS ; fi
if ! [ $N_JOBS -eq $N_JOBS ] 2>/dev/null ; then
    echo "ERROR: invalid number of jobs '$N_JOBS'."
    exit 1
fi

RELEASE=$3
if [ "x$RELEASE" = "x" ] ; then
    if [ $MODE = "limits" ] ; then
        RELEASE=$DEFAULT_RELEASE_LIMITS
    elif [ $MODE = "prod" ] ; then
        RELEASE=$DEFAULT_RELEASE_PROD
    elif [ $MODE = "prod17" ] ; then
        RELEASE=$DEFAULT_RELEASE_PROD17
    else
        RELEASE=$DEFAULT_RELEASE_ANA
    fi
fi
if [ -e $RELEASE ] ; then
echo "ERROR: Working area for $RELEASE already exists."
exit 1
fi

if [ $MODE = "limits" ] ; then
    export SCRAM_ARCH=slc6_amd64_gcc491
elif [ $MODE = "prod17" ] ; then
    export SCRAM_ARCH=slc6_amd64_gcc630
else
    export SCRAM_ARCH=slc6_amd64_gcc530
fi

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

if [ $MODE = "prod" ] ; then
    git cms-init

    # MET filters
    #git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate #outdated
    #git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter

    # MET corrections
    #git cms-merge-topic cms-met:METRecipe_8020
    git cms-merge-topic cms-met:METRecipe_8020_for80Xintegration
fi

if [ $MODE = "prod17" ] ; then
    git cms-init

    # Electron MVA identification
    git cms-merge-topic guitargeek:ElectronID_MVA2017_940pre3
    scram b -j 8

    # Add the area containing the MVA weights (from cms-data, to appear in “external”).
    # Note: the “external” area appears after “scram build” is run at least once, as above
    cd $RELEASE/external

    cd slc6_amd64_gcc630/
    git clone https://github.com/lsoffi/RecoEgamma-PhotonIdentification.git data/RecoEgamma/PhotonIdentification/data
    cd data/RecoEgamma/PhotonIdentification/data
    git checkout CMSSW_9_4_0_pre3_TnP
    cd $RELEASE/external
    cd slc6_amd64_gcc630/
    git clone https://github.com/lsoffi/RecoEgamma-ElectronIdentification.git data/RecoEgamma/ElectronIdentification/data
    cd data/RecoEgamma/ElectronIdentification/data
    git checkout CMSSW_9_4_0_pre3_TnP
    # Go back to the src/
    cd $RELEASE/src

fi

if [ $MODE = "limits" ] ; then
    # Combine tool
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git checkout v6.3.0
    cd ../..

    # CombineHarvester package
    git clone https://github.com/cms-analysis/CombineHarvester.git
    # HH stat tools
    git clone git@github.com:hh-italian-group/HHStatAnalysis.git
    cd HHStatAnalysis
    git checkout ttbb-it
    cd ..
fi

# SVfit packages
#git clone git@github.com:veelken/SVfit_standalone.git TauAnalysis/SVfitStandalone #notworking
git clone git@github.com:hh-italian-group/SVfit_standalone.git TauAnalysis/SVfitStandalone
cd TauAnalysis/SVfitStandalone
#git checkout HIG-16-006
git checkout hh_italian
cd ../..

# HHKinFit2 packages
git clone git@github.com:hh-italian-group/HHKinFit2.git HHKinFit2/HHKinFit2

# LeptonEfficiencies packages
git clone git@github.com:hh-italian-group/LeptonEff-interface.git HTT-utilities
#git clone git@github.com:CMS-HTT/LeptonEfficiencies.git HTT-utilities/LepEffInterface/data
git clone git@github.com:hh-italian-group/LeptonEfficiencies.git HTT-utilities/LepEffInterface/data

# Recoil Corrections
if [ $MODE = "prod" -o $MODE = "prod17" -o $MODE = "limits" ] ; then
    git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
fi

# hh-italian-group packages
git clone git@github.com:hh-italian-group/AnalysisTools.git
git clone git@github.com:hh-italian-group/h-tautau.git
git clone git@github.com:hh-italian-group/hh-bbtautau.git

if [ $MODE = "prod" ] ; then
    cd AnalysisTools
    git checkout prod_v3
    cd ..
    cd h-tautau
    git checkout prod_v3
    cd ..
    cd hh-bbtautau
    git checkout ana_v2
    cd ..
fi

if [ $MODE = "prod17" ] ; then
    cd AnalysisTools
    git checkout master
    cd ..
    cd h-tautau
    git checkout prod_v4_2017
    cd ..
    cd hh-bbtautau
    git checkout ana_v4
    cd ..
fi

if [ $MODE = "limits" ] ; then
    cd h-tautau
    git checkout sync
    cd ..
fi

if [ $MODE = "ana" ] ; then
    cd AnalysisTools
    git checkout master
    cd ..
    cd h-tautau
    git checkout ana_v3
    cd ..
    cd hh-bbtautau
    git checkout ana_v3
    cd ..
fi

# Prepare analysis working area
./AnalysisTools/Run/install.sh ../build AnalysisTools h-tautau hh-bbtautau

if [ $MODE = "prod" -o $MODE = "prod17" -o $MODE = "limits" ] ; then
    scram b -j$N_JOBS
fi
