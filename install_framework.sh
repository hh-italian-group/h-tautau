#!/usr/bin/env bash
# Install hh-italian-group framework.
# This file is part of https://github.com/hh-italian-group/h-tautau.

declare -A INSTALL_MODES
INSTALL_MODES=( ["prod"]="CMSSW_10_2_20 _amd64_gcc700" \
                ["ana"]="CMSSW_11_1_0_pre6 _amd64_gcc820" \
                ["ana_osx"]="bbtautau None")
DEFAULT_N_JOBS=4

function join_by { local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}"; }

if [ $# -lt 1 -o $# -gt 3 ] ; then
    echo "Usage: mode [n_jobs] [cmssw_release]"
    printf "\n\t%-20sinstallation mode. Supported modes: " "mode"
    join_by ", " "${!INSTALL_MODES[@]}"
    printf ".\n\t%-20sthe number of jobs to run simultaneously during the compilation." "n_jobs"
    printf " Default: $DEFAULT_N_JOBS.\n"
    printf "\t%-20sCMSSW release." "cmssw_release"
    printf " Defaults:\n"
    for mode in "${!INSTALL_MODES[@]}" ; do
        MODE_DESC=( ${INSTALL_MODES[$mode]} )
        printf "\t%-20s\t%-10s\t%s\n" "" "$mode" "${MODE_DESC[0]}"
    done
    exit 1
fi

MODE=$1
if [[ ! ${!INSTALL_MODES[@]} =~ $MODE ]] ; then
    echo "ERROR: unsupported installation mode '$MODE'."
    printf "Supported installation modes: "
    join_by ", " "${!INSTALL_MODES[@]}"
    printf ".\n"
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
    MODE_DESC=( ${INSTALL_MODES[$MODE]} )
    RELEASE=${MODE_DESC[0]}
    if [ "$MODE" != "ana_osx" ] ; then
        export SCRAM_ARCH="slc$(cat /etc/redhat-release | sed -E 's/[^67]*([67])\..*/\1/')${MODE_DESC[1]}"
    fi
fi
if [ -e $RELEASE ] ; then
    echo "ERROR: Working area for $RELEASE already exists."
    exit 1
fi

function run_cmd {
    "$@"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error while rinning '$@'"
        exit 1
    fi
}

if [ "$MODE" != "ana_osx" ] ; then
    run_cmd scramv1 project CMSSW $RELEASE
    cd $RELEASE/src
    run_cmd eval `scramv1 runtime -sh`
else
    mkdir -p "$RELEASE"
    cd "$RELEASE"
fi

if [ $MODE = "prod" ] ; then
    run_cmd git cms-init
    run_cmd git cms-addpkg RecoMET/METFilters
    #run_cmd git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X
    run_cmd git cms-merge-topic cms-egamma:EgammaPostRecoTools

	# Update PileupJetID
	run_cmd git cms-addpkg  RecoJets/JetProducers
	run_cmd git clone -b 94X_weights_DYJets_inc_v2 git@github.com:cms-jet/PUjetID.git PUJetIDweights/
	run_cmd cp PUJetIDweights/weights/pileupJetId_{94,102}X_Eta* $CMSSW_BASE/src/RecoJets/JetProducers/data/
	run_cmd rm -rf PUJetIDweights/
	run_cmd git cms-merge-topic alefisico:PUID_102X
fi

# new SVfit packages
run_cmd git clone -b hh-italian git@github.com:hh-italian-group/ClassicSVfit.git TauAnalysis/ClassicSVfit
run_cmd git clone git@github.com:hh-italian-group/SVfitTF.git TauAnalysis/SVfitTF

# HHKinFit2 packages
run_cmd git clone git@github.com:hh-italian-group/HHKinFit2.git HHKinFit2/HHKinFit2

# LeptonEfficiencies packages
run_cmd git clone git@github.com:hh-italian-group/LeptonEff-interface.git HTT-utilities
run_cmd git clone git@github.com:hh-italian-group/LeptonEfficiencies.git HTT-utilities/LepEffInterface/data

# Tau ID and Trigger SFs
run_cmd git clone git@github.com:cms-tau-pog/TauIDSFs.git TauPOG/TauIDSFs
run_cmd git clone -b run2_SFs git@github.com:cms-tau-pog/TauTriggerSFs.git TauAnalysisTools/TauTriggerSFs

run_cmd git clone git@github.com:hh-italian-group/HHbtag.git HHTools/HHbtag

# Recoil Corrections
#run_cmd git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections

# Install analysis packages
declare -A ANA_PACKAGES
ANA_PACKAGES=( ["AnalysisTools"]="prod:master ana:master ana_osx:master" \
               ["h-tautau"]="prod:master ana:master ana_osx:master" \
               ["hh-bbtautau"]="prod:None ana:master ana_osx:master" )
GITHUB_USER=$(git config user.github)

for pkg in "${!ANA_PACKAGES[@]}" ; do
    pkg_descs="${ANA_PACKAGES[$pkg]}"
    branch="master"
    for desc in $pkg_descs ; do
        if [ "${desc%%:*}" = "$MODE" ] ; then
            branch=${desc##*:}
            break
        fi
    done

    if [ "$branch" = "None" ] ; then continue; fi

    run_cmd git clone git@github.com:hh-italian-group/${pkg}.git
    cd "$pkg"
    if [ "$branch" != "master" ] ; then
        run_cmd git checkout -b $branch origin/$branch
    fi
    git ls-remote git@github.com:$GITHUB_USER/${pkg}.git &> /dev/null
    RESULT=$?
    if [ $RESULT -eq 0 ] ; then
        git remote add $GITHUB_USER git@github.com:$GITHUB_USER/${pkg}.git
        git fetch $GITHUB_USER
    fi
    cd ..
done

if [ "$MODE" = "ana_osx" ] ; then
    BUILD_PATH=build
else
    BUILD_PATH=../build
fi

# Prepare analysis working area
run_cmd ./AnalysisTools/Run/install.sh "$BUILD_PATH" "${!ANA_PACKAGES[@]}"

if ! [[ $MODE =~ ana.* ]] ; then
    run_cmd scram b -j$N_JOBS
fi

echo "Framework has been successfully installed."
