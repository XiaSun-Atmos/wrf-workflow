#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Date: Jul, 2023


module load hpss

ulimit -s unlimited
ulimit -u 9600

# Archive locations
SAVE_BDY=/ESRL/BMC/wrfruc/5year/Xia.Sun/wfip3_cases
# SAVE_WRFOUT=/mnt/lfs3/projects/nrtrr/wrfout/hrrr
SAVE_UPP=/ESRL/BMC/wrfruc/5year/Xia.Sun/wfip3_cases
DISK_BDY=/home/Xia.Sun/lfs5_rtwbl/wfip3_cases/data/IC_BC

DATAHOME=${RUN_DIR}/run_$YYYY$MM$DD$HH

# Check to make sure that the DATAHOME exists
if [ ! ${DATAHOME} ]; then
  ${ECHO} "ERROR: DATAHOME, \$DATAHOME, is not defined"
  exit 1
fi

# Archive analysis (initial condition) and boundary condition
echo 'Save IC for '${YYYY}${MM}${DD}${HH}
if [ -r ${DATAHOME}/wrfinput_d01 ] && [ -r ${DATAHOME}/wrfbdy_d01 ]; then

  htar -cvf ${SAVE_BDY}/wrfinput_${YYYY}${MM}${DD}${HH}_exp2.tar ${DATAHOME}/wrfinput_d01 ${DATAHOME}/wrfinput_d02 ${DATAHOME}/wrfbdy_d01 ${DATAHOME}/namelist.input
fi 


# Archive WRF output
# htar -cvf ${SAVE_BDY}/wrfout_d01_${YYYY}${MM}${DD}${HH}.exp1.tar ${DATAHOME}/wrfout_d01*
# htar -cvf ${SAVE_BDY}/wrfout_d02_${YYYY}${MM}${DD}${HH}.exp1.tar ${DATAHOME}/wrfout_d02*

# Archive GRIB2 output
echo 'Save GRIB2 files for '${YYYY}${MM}${DD}${HH}
htar -cvf ${SAVE_UPP}/${YYYY}${MM}${DD}${HH}.grib2.exp2.tar ${RUN_DIR}/$YYYY$MM$DD$HH/
