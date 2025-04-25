#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Date: Feb, 2023

limit stacksize unlimited

# clean up WPS dir
rm -r ${RUN_DIR}/WPS_$YYYY$MM$DD$HH/FILE:202*
rm -r ${RUN_DIR}/WPS_$YYYY$MM$DD$HH/met_em.d0*.nc
rm -rf ${RUN_DIR}/WPS_$YYYY$MM$DD$HH

DATAHOME=${RUN_DIR}/run_$YYYY$MM$DD$HH
# move IC/BC
DISK_BDY=/home/Xia.Sun/lfs5_rtwbl/wfip3_cases/data/IC_BC
mkdir ${DISK_BDY}/$YYYY$MM$DD$HH
cp -r ${DATAHOME}/wrfinput_d0* ${DISK_BDY}/$YYYY$MM$DD$HH/
cp -r ${DATAHOME}/wrfbdy_d01 ${DISK_BDY}/$YYYY$MM$DD$HH/
cp -r ${DATAHOME}/namelist.input ${DISK_BDY}/$YYYY$MM$DD$HH/
cp -r ${DATAHOME}/wrfqnain* ${DISK_BDY}/$YYYY$MM$DD$HH/
cp -r ${DATAHOME}/wrflowinp* ${DISK_BDY}/$YYYY$MM$DD$HH/

# move WRF outputs to its designated location
mkdir ${OUT_DIR}/wrf_out/
mkdir ${OUT_DIR}/wrf_out/$YYYY$MM$DD$HH
mv ${RUN_DIR}/run_$YYYY$MM$DD$HH/wrfout* ${OUT_DIR}/wrf_out/$YYYY$MM$DD$HH/
# rm -rf ${RUN_DIR}/run_$YYYY$MM$DD$HH

# move UPP outputs to its designated location
mv ${RUN_DIR}/$YYYY$MM$DD$HH ${OUT_DIR}/upp_out/

# clean up GFS IC dir
# rm -rf ${GFS_DIR}/${YYYY}${MM}${DD}${HH}00
exit