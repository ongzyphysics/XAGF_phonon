#!/bin/bash

CURRDIR=$(pwd)
AGFDIR=${CURRDIR} # working directory containing XAGF code
INPDIR=${CURRDIR}/input_files_cnt  # directory containing input files for AGF calculations
OUTDIR=${CURRDIR}/output_files_cnt # directory containing output files from AGF calculations

cd ${AGFDIR}
matlab -nodisplay -nodesktop -r "ExtendedAtomisticGreensFunctionTransmission('${INPDIR}','${OUTDIR}'); exit" 
cd ${CURRDIR}

exit 1

# ===================================================================================

