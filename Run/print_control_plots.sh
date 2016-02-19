#!/bin/bash
# Print all control plots to put to the AN.
# This file is part of https://github.com/hh-italian-group/h-tautau.

if [ $# -ne 2 ] ; then
    echo "Usage: input_file output_path"
    exit
fi

INPUT_FILE=$1
OUTPUT_PATH=$2

SOURCE_CFG="Analysis/config/sources.cfg"
CHANNEL_NAME="tauTau"
SIGNAL_NAME="ggHhh300"
IS_BLIND=false
DRAW_RATIO=false
DRAW_BKG_ERRORS=true
USE_LOG_Y=false

SEL_2JET="2jets/NoCuts/OS_Isolated/Central"
SEL_0TAG_NO_MASS_CUT="2jets0btag/NoCuts/OS_Isolated/Central"
SEL_1TAG_NO_MASS_CUT="2jets1btag/NoCuts/OS_Isolated/Central"
SEL_2TAG_NO_MASS_CUT="2jets2btag/NoCuts/OS_Isolated/Central"
SEL_0TAG_MASS_CUT="2jets0btag/MassWindow/OS_Isolated/Central"
SEL_1TAG_MASS_CUT="2jets1btag/MassWindow/OS_Isolated/Central"
SEL_2TAG_MASS_CUT="2jets2btag/MassWindow/OS_Isolated/Central"

if [ ! -f "$INPUT_FILE" ] ; then
    echo "ERROR: input file '$INPUT_FILE' not found."
    exit
fi

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' not found."
    exit
fi

MAKE_PATH="./RunTools/make_withFactory.sh"
CMD_NAME="Print_Stack"
EXE_NAME="$OUTPUT_PATH/$CMD_NAME"

$MAKE_PATH $OUTPUT_PATH $CMD_NAME $CMD_NAME
RESULT=$?
if [ $RESULT -ne 0 ] ; then
    echo "Compilation of $EXE_NAME failed with an error code $RESULT."
    exit
fi

function print_plot {
    local PLOT_NAME=$1
    local SELECTION=$2
    local OUTPUT_DIR=$3
    local OUTPUT_NAME=$4
    echo "Printing $PLOT_NAME in $SELECTION ..."
    $EXE_NAME $SOURCE_CFG $INPUT_FILE $OUTPUT_PATH $CHANNEL_NAME $SELECTION $PLOT_NAME $SIGNAL_NAME \
              "@$IS_BLIND" "@$DRAW_RATIO" "@$DRAW_BKG_ERRORS"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        exit
    fi
    local SELECTION_SUFFIX=$( echo $SELECTION | sed 's/\//_/g' )
    local RESULT_FILE_PATH="$OUTPUT_PATH/${CHANNEL_NAME}_${PLOT_NAME}_${SELECTION_SUFFIX}.pdf"
    if [ ! -f "$RESULT_FILE_PATH" ] ; then
        echo "ERROR: Output file '$RESULT_FILE_PATH' not found."
        exit
    fi
    local OUTPUT_DIR_PATH="$OUTPUT_PATH/$OUTPUT_DIR"
    mkdir -p $OUTPUT_DIR_PATH
    local OUTPUT_FILE_PATH="$OUTPUT_DIR_PATH/$OUTPUT_NAME.pdf"
    mv $RESULT_FILE_PATH $OUTPUT_FILE_PATH
}

print_plot "MET" $SEL_2JET "ControlPlots/tautau" "inclusive_MET"
print_plot "csv_b1" $SEL_2JET "ControlPlots/tautau" "inclusive_csv_b1"
print_plot "csv_b2" $SEL_2JET "ControlPlots/tautau" "inclusive_csv_b2"
print_plot "eta_1" $SEL_2JET "ControlPlots/tautau" "inclusive_eta_1"
print_plot "eta_2" $SEL_2JET "ControlPlots/tautau" "inclusive_eta_2"
print_plot "eta_b1" $SEL_2JET "ControlPlots/tautau" "inclusive_eta_b1"
print_plot "eta_b2" $SEL_2JET "ControlPlots/tautau" "inclusive_eta_b2"
print_plot "pt_1" $SEL_2JET "ControlPlots/tautau" "inclusive_pt_1"
print_plot "pt_2" $SEL_2JET "ControlPlots/tautau" "inclusive_pt_2"
print_plot "pt_b1" $SEL_2JET "ControlPlots/tautau" "inclusive_pt_b1"
print_plot "pt_b2" $SEL_2JET "ControlPlots/tautau" "inclusive_pt_b2"

#print_plot "m_ttbb_kinfit" $SEL_0TAG_MASS_CUT "MassPlots/tautau" "0tag_m_H_kinfit_masscut"
print_plot "m_ttbb_kinfit" $SEL_0TAG_NO_MASS_CUT "MassPlots/tautau" "0tag_m_H_kinfit_no_mass_cut"
print_plot "m_ttbb" $SEL_0TAG_NO_MASS_CUT "MassPlots/tautau" "0tag_m_H_no_kinfit_no_mass_cut"
print_plot "m_bb" $SEL_0TAG_NO_MASS_CUT "MassPlots/tautau" "0tag_m_bb"
print_plot "m_sv" $SEL_0TAG_NO_MASS_CUT "MassPlots/tautau" "0tag_m_sv"

#print_plot "m_ttbb_kinfit" $SEL_1TAG_MASS_CUT "MassPlots/tautau" "1tag_m_H_kinfit_masscut"
print_plot "m_ttbb_kinfit" $SEL_1TAG_NO_MASS_CUT "MassPlots/tautau" "1tag_m_H_kinfit_no_mass_cut"
print_plot "m_ttbb" $SEL_1TAG_NO_MASS_CUT "MassPlots/tautau" "1tag_m_H_no_kinfit_no_mass_cut"
print_plot "m_bb" $SEL_1TAG_NO_MASS_CUT "MassPlots/tautau" "1tag_m_bb"
print_plot "m_sv" $SEL_1TAG_NO_MASS_CUT "MassPlots/tautau" "1tag_m_sv"

#print_plot "m_ttbb_kinfit" $SEL_2TAG_MASS_CUT "MassPlots/tautau" "2tag_m_H_kinfit_masscut"
print_plot "m_ttbb_kinfit" $SEL_2TAG_NO_MASS_CUT "MassPlots/tautau" "2tag_m_H_kinfit_no_mass_cut"
print_plot "m_ttbb" $SEL_2TAG_NO_MASS_CUT "MassPlots/tautau" "2tag_m_H_no_kinfit_no_mass_cut"
print_plot "m_bb" $SEL_2TAG_NO_MASS_CUT "MassPlots/tautau" "2tag_m_bb"
print_plot "m_sv" $SEL_2TAG_NO_MASS_CUT "MassPlots/tautau" "2tag_m_sv"


rm -f $EXE_NAME
