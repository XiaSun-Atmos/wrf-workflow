#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# History:
# Transitioned to Rocky8 from CentOS on Mar 13, 2024 with updated modules
# Established Date: Feb, 2023
# 

# load modules on Jet needed for WPS exes

module purge
module load intel/2022.1.2
module load impi/2022.1.2
module load netcdf/4.7.0
module load szip/2.1  hdf5/1.10.5
module load netcdf/4.7.0
module load pnetcdf/1.10.0
module load jasper/4.2.0
# export MPIRUN="mpirun"
export MPIRUN="srun"
limit stacksize unlimited
ulimit -s unlimited
ulimit -u 9600
# set up start_date and end_date for namelist.input
rm -rf ${RUN_DIR}/WPS_${YYYY}${MM}${DD}${HH}
cp -r ${WPS_DIR}/WPS_Para ${RUN_DIR}/WPS_${YYYY}${MM}${DD}${HH}
cd ${RUN_DIR}/WPS_${YYYY}${MM}${DD}${HH}


cp -r ${NML_DIR}/namelist.wps.wfip3 namelist.wps
sed -i "s/STARTYY/$YYYY/g" namelist.wps
sed -i "s/STARTMM/$MM/g" namelist.wps
sed -i "s/STARTDD/$DD/g" namelist.wps
sed -i "s/STARTHH/$HH/g" namelist.wps
end_date_str=$(date -u -d "${YYYY}-${MM}-${DD} ${HH}:00:00 ${fcsthrs} hours" "+%Y-%m-%d_%H:00:00")
sed -i "s/ENDDATE/${end_date_str}/g" namelist.wps


# ./geogrid.exe >& log.geogrid_${YYYY}-${MM}-${DD}-${HH}

if [ -e "geo_em.d02.nc" ]; then
  echo "Success: geo_em.d02* file generated"
else
  echo "Error: geo_em.d02* file not generated"
fi

./link_grib.csh ${GFS_DIR}/${YYYY}${MM}${DD}${HH}/*
# ln -sf ungrib/Variable_Tables/Vtable.GFS Vtable
ln -sf Vtable.RAP Vtable
./ungrib.exe >& ungrib.exe.log_${YYYY}-${MM}-${DD}-${HH}
echo "ungrib.exe finished"

# link OSTIA intermediate files
ln -s /lfs6/BMC/rtwbl/Xia.Sun/wfip3/data/ostia/*/SST* .

# Now run metgrid
srun ./metgrid.exe >& log.metgrid_${YYYY}-${MM}-${DD}-${HH}

# $MPIRUN metgrid.exe >& metgrid.exe.log_${YYYY}-${MM}-${DD}-${HH}

if [ -e "met_em.d02.${YYYY}-${MM}-${DD}_00:00:00.nc" ]; then
  echo "Success: met_em.d02* file generated"
else
  echo "met_em.d02.${YYYY}-${MM}-${DD}_00:00:00.nc"
  echo "Error: met_em.d02* file not generated"
fi
