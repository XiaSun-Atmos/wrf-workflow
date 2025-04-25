#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# History:
# Transitioned to Rocky8 from CentOS on Mar 13, 2024 with updated modules
# Established Date: Feb, 2023

# load modules on Jet needed for WRF exes
#module purge
#module load intel/2022.1.2
#module load impi/2022.1.2
#module load netcdf/4.7.0
#module load szip/2.1  hdf5/1.10.5
#module load netcdf/4.7.0
#module load pnetcdf/1.10.0
#module load jasper/4.2.0

# Set IMPI I/O performance variables
#export I_MPI_EXTRA_FILESYSTEM=on
#export I_MPI_EXTRA_FILESYSTEM_LIST=lustre:panfs
#export OMP_NUM_THREADS=4
ulimit -s unlimited
module purge
module load gnu
module load intel/2023.2.0
module load impi/2023.2.0
module load pnetcdf/1.12.3
module load netcdf
module load hdf5/1.14.3
export PNETCDF=/apps/pnetcdf/1.12.3/intel_2023.2.0-impi/
# running wrf.exe
cd ${RUN_DIR}/run_$YYYY$MM$DD$HH

cp -L ${NML_DIR}/namelist.input.wfip3.rrfs namelist.input

sed -i "s/STARTYY/$YYYY/g" namelist.input
sed -i "s/STARTMM/$MM/g" namelist.input
sed -i "s/STARTDD/$DD/g" namelist.input
sed -i "s/STARTHH/$HH/g" namelist.input


end_date_str=$(date -u -d "${YYYY}-${MM}-${DD} ${HH}:00:00 ${fcsthrs} hours" "+%Y-%m-%d_%H:00:00")
end_year=$(echo "$end_date_str" | cut -d '-' -f 1)
end_month=$(echo "$end_date_str" | cut -d '-' -f 2)
end_day=$(echo "$end_date_str" | cut -d '-' -f 3 | cut -d '_' -f 1)
end_hour=$(echo "$end_date_str" | cut -d '_' -f 2 | cut -d ':' -f 1)

sed -i "s/ENDYY/${end_year}/g" namelist.input
sed -i "s/ENDMM/${end_month}/g" namelist.input
sed -i "s/ENDDD/${end_day}/g" namelist.input
sed -i "s/ENDHH/${end_hour}/g" namelist.input

echo ${fcsthrs}
# export MPIRUN="mpirun"
export MPIRUN="srun"
export WRF=${RUN_DIR}/run_$YYYY$MM$DD$HH/wrf.exe
$MPIRUN $WRF >& wrfexe.log_${YYYY}-${MM}-${DD}-${HH}
