#!/bin/bash

CURRDIR=$(pwd)
AGFDIR=${CURRDIR}
INPDIR=${CURRDIR}/input_files_cnt
OUTDIR=${CURRDIR}/output_files_cnt
# AGFDIR=/home/ongzy/IHPC_Documents/ihpc_research/2017-12-12_FullSMatrix/SourceCode
# INPDIR=/home/ongzy/IHPC_Documents/ihpc_research/2018-08-29_CNTJunction/input_files
# OUTDIR=/home/ongzy/IHPC_Documents/ihpc_research/2018-08-29_CNTJunction/output_files

cd ${AGFDIR}
matlab -nodisplay -nodesktop -r "ExtendedAtomisticGreensFunctionTransmission('${INPDIR}','${OUTDIR}'); exit" 
cd ${CURRDIR}

# ===================================================================================
exit 1


# ===================================================================================

cp ${CURRDIR}/IMJ_4_pairs/input_files/{Left,Center,Right}_*.agf ${CURRDIR}/input_files_lowfreq/

cd ${AGFDIR}
matlab -nodisplay -nodesktop -r "ExtendedAtomisticGreensFunctionTransmission('${INPDIR}','${OUTDIR}'); exit" 
cd ${CURRDIR}

cp ${CURRDIR}/output_files_lowfreq/*.dat ${CURRDIR}/IMJ_4_pairs/output_files_lowfreq/

# ===================================================================================
exit 1

# ===================================================================================

cp ${CURRDIR}/IMJ_8_pairs/input_files/{Left,Center,Right}_*.agf ${CURRDIR}/input_files_lowfreq/

cd ${AGFDIR}
matlab -nodisplay -nodesktop -r "ExtendedAtomisticGreensFunctionTransmission('${INPDIR}','${OUTDIR}'); exit" 
cd ${CURRDIR}

cp ${CURRDIR}/output_files_lowfreq/*.dat ${CURRDIR}/IMJ_8_pairs/output_files_lowfreq/
