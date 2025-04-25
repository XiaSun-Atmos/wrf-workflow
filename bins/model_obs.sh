#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Date: Dec, 2024
# Purpose: Plotting time height comparisons and surface met comparisons 
# Reference: https://www.psl.noaa.gov/renewable_energy/wfip3/modelobs/



limit stacksize unlimited

module load cdo
module load nco

column_d1_dir=${OUT_DIR}/column_out/${YYYY}${MM}${DD}${HH}/d01/output
column_d2_dir=${OUT_DIR}/column_out/${YYYY}${MM}${DD}${HH}/d02/output
model_dir=${OUT_DIR}/column_out/wfip3_experimental
d1_dir=${OUT_DIR}/column_out/wfip3_experimental/3km/${YYYY}
d2_dir=${OUT_DIR}/column_out/wfip3_experimental/1km/${YYYY}
psl_dir=${OUT_DIR}/column_out


hrrr_dir=${OUT_DIR}/column_out/rapV5_hrrrV4/hrrr/${YYYY}

export OUT_DIR=${OUT_DIR}
export BIN_DIR=${BIN_DIR}
export PSL_DIR=${psl_dir}
export YYYY=${YYYY}
export MM=${MM}
export DD=${DD}
export HH=${HH}
# Check if the directory exists


if [ ! -d "$d1_dir" ]; then
    # Create the directory
    mkdir -p "$d1_dir"
    echo "Directory '$d1_dir' created."
else
    echo "Directory '$d1_dir' existing"
fi

if [ ! -d "$d2_dir" ]; then
    # Create the directory
    mkdir -p "$d2_dir"
    echo "Directory '$d2_dir' created."
else
    echo "Directory '$d2_dir' existing"
fi

if [ ! -d "$psl_dir/Measurement_Data/Stations" ]; then
    # Create the directory
    mkdir -p "${psl_dir}/Measurement_Data/Stations"
    echo "Directory '$psl_dir' created."
    cp -r /lfs5/BMC/rtwbl/Xia.Sun/wfip3_cases/data/obs/Measurement_Data/Stations/* ${psl_dir}/Measurement_Data/Stations/
else
    echo "Directory '$psl_dir' existing"
fi


# concatenate hourly profiler.nc files
first_file=$(ls ${column_d1_dir}/profiler_sitesF.ncep_wfip3km* | head -n 1 | xargs -n 1 basename)
doy=$(echo "${first_file}" | grep -oP '(?<=wfip3km\.)[^.]+')


hrrr_file="profiler_sitesF.ncep_hrrr.${doy}.${HH}00.nc"

if [ ! -f "$hrrr_dir/${hrrr_file}" ]; then
    # Create the directory
    mkdir -p "$hrrr_dir"
    echo "Directory '$hrrr_dir' created."
    cp -r /lfs5/BMC/rtwbl/Xia.Sun/wfip3_cases/outputs/column_out/rapV5_hrrrV4/hrrr/${YYYY}/${hrrr_file} ${hrrr_dir}/
else
    echo "Directory '$hrrr_dir' existing"
fi

ncrcat -O -h ${column_d1_dir}/profiler_sitesF.ncep_wfip3km.*.00.nc ${d1_dir}/profiler_sitesF.ncep_wfip3.${doy}.${HH}00.nc
# if running 2 domains, then uncomment below; Xia Sun Dec 11 2024
ncrcat -O -h ${column_d2_dir}/profiler_sitesF.ncep_wfip1km.*.00.nc ${d2_dir}/profiler_sitesF.ncep_wfip3.${doy}.${HH}00.nc



module use /contrib/miniconda3/modulefiles
module load miniconda3
conda activate /lfs5/BMC/rtwbl/Xia.Sun/envs/py31013

python ${BIN_DIR}/model_obs/Main_Driver/drivers2run/rwp_driver.py 

python ${BIN_DIR}/model_obs/Main_Driver/drivers2run/pslmet_driver.py
