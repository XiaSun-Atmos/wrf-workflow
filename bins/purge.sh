#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Purge GRIB and wrfout files after two days
# Date: Jul, 2023


limit stacksize unlimited

# clean up WRF outputs after two days on disk
rm -rf ${OUT_DIR}/wrf_out/$YYYY$MM$DD$HH

echo "cleaning WRF outputs"
# clean up UPP outputs after two days on disk
rm -rf ${OUT_DIR}/upp_out/$YYYY$MM$DD$HH

echo "cleaning UPP outputs"
# clean up log files and stdout files on disk
rm -r ${STDOUT}/*/$MM$DD$HH*

echo "cleaning logfiles and stdout files"
#clean up log files
rm -r ${LOG}/workflow_$YYYY$MM$DD$HH*
rm -r ${LOG}/python_$YYYY$MM$DD$HH*

#clean up IC
rm -r ${GFS_DIR}/$YYYY$MM$DD$HH
echo "cleaning IC files"

#clean up python graph zip files
rm -rf ${OUT_DIR}/graphs/$YYYY$MM$DD$HH
echo "cleaning python graph zip files"