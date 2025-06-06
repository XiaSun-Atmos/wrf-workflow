#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Date: Sep, 2024
# Purpose: Subset WRF outputs with specified lat and lon ranges

limit stacksize unlimited

module load cdo

WRFDIR=${OUT_DIR}/wrf_out/$YYYY$MM$DD$HH/
SUBDIR=${OUT_DIR}/wrf_subset/$YYYY$MM$DD$HH/

mkdir ${OUT_DIR}/wrf_subset/$YYYY$MM$DD$HH

InitTime="${YYYY}-${MM}-${DD} ${HH}:00:00"
StartTime=$(date -d "$InitTime 6 hours" +"%Y-%m-%d %H:00:00")
EndTime=$(date -d "$InitTime $fcsthrs hours" +"%Y-%m-%d_%H:00:00")

interval_minutes=15

total_intervals=$((${fcsthrs} * 60 / interval_minutes))

export WRFDIR
export SUBDIR
export StartTime
export interval_minutes

wrf_subset (){
    # Add 15 * i minutes to the original date
     echo ${1}

    new_date=$(date -d "$StartTime $((${1} * interval_minutes)) minutes" +"%Y-%m-%d_%H_%M_%S")
    
    # Output the new date
    echo $new_date

    # process wrfouts
 #cdo sellonlatbox,-72.7701,-69.3134,40.0273,42.098 ${WRFDIR}/wrfout_d01_${new_date} ${SUBDIR}/wrfout_d01_${new_date}_subset
 #cdo sellonlatbox,-72.7701,-69.3134,40.0273,42.098 ${WRFDIR}/wrfout_d02_${new_date} ${SUBDIR}/wrfout_d02_${new_date}_subset
  # extend 50 km each direction from the original area
 cdo sellonlatbox,-73.2199,-68.8636,39.5773,42.5480 ${WRFDIR}/wrfout_d01_${new_date} ${SUBDIR}/wrfout_d01_${new_date}_subset
 cdo sellonlatbox,-73.2199,-68.8636,39.5773,42.5480 ${WRFDIR}/wrfout_d02_${new_date} ${SUBDIR}/wrfout_d02_${new_date}_subset
}
export -f wrf_subset
parallel -u wrf_subset ::: {0..168}
