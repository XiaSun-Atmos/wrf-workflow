#!/bin/bash
# Mar, 2024: Transitioned to Rocky8 from CentOS
# Jul, 2023: Xia Sun added sub-hourly capbility
set -x

module use /lfs5/BMC/rtwbl/Xia.Sun/upp_exes/UPP/modulefiles
module load jet
module load stack-intel/2021.5.0 #Xia 2025-01-09
module load sfcio                #Xia 2025-01-09
module load w3nco


#--------------------------------------------------------
# Updates:
#
# August 2005: Hui-Ya Chuang, NCEP: This script uses 
# NCEP's Unipost to post processes WRF native model 
# output, and uses copygb to horizontally interpolate posted 
# output from native A-E to a regular projection grid. 
#
# July 2006: Meral Demirtas, NCAR/DTC: Added new "copygb" 
# options and revised some parts for clarity. 
#
# April 2015: Modified to run NMM-B/NEMS, KRF(DTC)
#
# January 2019: Modified to remove NMMB/NMM and add FV3, TH (DTC)
# 
# October 2019: Modified for new unified build system;
#               ndate.exe and copygb.exe have been removed
#
# October 2020: Modified to remove WRF and grib1; Add FV3LAM
#               Updates for cmake build, Change exec name             
#
# May 2022: Modified to remove binarynemsiompiio; 
#           Added netcdfpara; Removed netcdf;
#           Changed UPP directory name and path             
#
# June 2022: Add 2D decomposition capabilities, TH (DTC)
#--------------------------------------------------------
#
# This script runs the stand-alone community version of UPP
#
#--------------------------------------------------------

#----------------------------------------------------------------------------------
#--- USER EDIT DESCIPTIONS --------------------------------------------------------
# See UPP User's Guide for more information
# https://upp.readthedocs.io/en/latest/
#----------------------------------------------------------------------------------
# TOP_DIR       : Top level directory for building and running UPP
# DOMAINPATH    : Working directory for this run.
# UPP_HOME      : Location of the UPP directory
# POSTEXEC      : Location of the UPP executable
# modelDataPath : Location of the model output data files to be post-processed
# txtCntrlFile  : Name and location of the flat text file that lists desired fields for output
#                 GFS: postxconfig-NT-GFS-F00.txt (0-hour lead) and postxconfig-NT-GFS.txt (all other
#                 leads)
#                 LAM (Limited Area Model): postxconfig-NT-fv3lam.txt
# model         : What model is used? GFS or LAM (Limited Area Model)
# inFormat      : Format of the model data 
#                 GFS - "netcdfpara"
#                 LAM - "netcdfpara"
# outFormat     : Format of output from UPP 
#                 grib2 
# startdate     : Forecast start date (YYYYMMDDHH)
# fhr           : First forecast hour to be post-processed
# lastfhr       : Last forecast hour to be post-processed
# incrementhr   : Increment (in hours) between forecast files
#                 * Do not set to 0 or the script will loop continuously *
# RUN_COMMAND   : System run command for serial or parallel runs, examples below.
#
# numx          : The number of subdomains in the x-direction for 2D decomposition.
#                 (default is 1)
#
#----------------------------------------------------------------------------------
#--- BEGIN USER EDIT HERE ---------------------------------------------------------
#----------------------------------------------------------------------------------
DATE=/bin/date
# Set relevant paths and data information
# This script assumes you created a directory $DOMAINPATH/postprd
# as recommended in the users guide where UPP will output.
export TOP_DIR=/mnt/lfs5/BMC/rtwbl/Xia.Sun/wfip3_cases/run_dir
export DOMAINPATH=${TOP_DIR}
# export UPP_HOME=${TOP_DIR}/upp/UPP
export UPP_HOME=${UPP_DIR}
export POSTEXEC=${UPP_HOME}/tests/install/bin
export modelDataPath=${RUN_DIR}/run_$YYYY$MM$DD$HH
#export txtCntrlFile=${DOMAINPATH}/postprd/parm/postxconfig-NT-WRF.txt
export txtCntrlFile=${NML_DIR}/postxconfig-NT-hrrr-highres.txt
# export txtCntrlFile=${DOMAINPATH}/postprd/parm/postxconfig-NT-wfip3.txt
# Specify model ("GFS" or "LAM" in upper case)
export model="RAPR"
# Set input format from model and ouput format from UPP
export inFormat="netcdfpara"
export outFormat="grib2"

# Set date/time information
export startdate=${YYYY}${MM}${DD}${HH}
export fhr=${fcsthr}
export lastfhr=${fcsthr}
export incrementhr=01
export startmin=${fcstmin}
export incrementmin=15
export lastmin=${fcstmin}
# Set run command: 
mkdir ${DOMAINPATH}/postprd_${startdate}_F${fhr}_${startmin}
cd ${DOMAINPATH}/postprd_${startdate}_F${fhr}_${startmin} && mkdir parm && cd parm
cp -r ${UPP_HOME}/parm/* ${DOMAINPATH}/postprd_${startdate}_F${fhr}_${startmin}/parm
rm -r  postxconfig-NT.txt
# Single processor  command example
#export RUN_COMMAND="${POSTEXEC}/upp.x "

# Parallel command examples:
#export RUN_COMMAND="mpirun -np 4 ${POSTEXEC}/upp.x "
export RUN_COMMAND="srun ${POSTEXEC}/upp.x "
#export RUN_COMMAND="mpirun.lsf ${POSTEXEC}/upp.x "
#export RUN_COMMAND="mpiexec_mpt ${POSTEXEC}/upp.x "

# The number of subdomains in the x-direction for 2d decomposition
export numx=2

# Shouldn't need to edit these.
# tmmark is an variable used as the file extention of the output
# filename .GrbF is used if this variable is not set
# COMSP is a variable used as the initial string of the output filename
export tmmark=tm00
export MP_SHARED_MEMORY=yes
export MP_LABELIO=yes

#----------------------------------------------------------------------
#--- END USER EDIT ----------------------------------------------------
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Shouldn't need to edit below unless something goes wrong or debugging
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Do some checks for directory/executable existence, user input, etc.
#----------------------------------------------------------------------
if [ ! -d ${POSTEXEC} ]; then
  echo "ERROR: POSTEXEC, '${POSTEXEC}', does not exist"
  exit 1
fi

if [ ! -x ${POSTEXEC}/upp.x ]; then
  echo "ERROR: upp.x, '${POSTEXEC}/upp.x', does not exist or is not executable."
  exit 1
fi

# Set tag based on user defined model (GFS or LAM)
if [ $model = "GFS" ]; then
   export tag=GFS
elif [ $model = "LAM" ]; then
   export tag=FV3R
else
   export tag=RAPR
fi

if [ ${model} == "GFS" ]; then
   if [ ${inFormat} == "netcdfpara" ]; then
      echo "Check: You are using 'model' 'inFormat'!"
   else
      echo "ERROR: 'inFormat' must be 'netcdfpara' for GFS model output. Exiting... "
      exit 1
   fi
elif [ ${model} == "LAM" ]; then
   if [ ${inFormat} == "netcdfpara" ]; then
      echo "Check: You are using 'model' 'inFormat'!"
   else
      echo "ERROR: 'inFormat' must be 'netcdfpara' for LAM  model output. Exiting... "
      exit 1
   fi
fi

if [[ ${outFormat} == "grib2" ]]; then
   if [ ! -e ${txtCntrlFile} ]; then
      echo "ERROR: 'txtCntrlFile' not found in '${txtCntrlFile}'.  Exiting... "
      exit 1
   fi
else
   echo "${outFormat} is not supported. Edit script to choose 'grib2' for 'outFormat'. Exiting... "
fi
 
if [ ! -d ${DOMAINPATH}/postprd_${startdate}_F${fhr}_${startmin} ]; then
  echo "ERROR: DOMAINPATH/postprd, '${DOMAINPATH}/postprd', does not exist. Exiting..."
  exit 1
fi

if [ ${incrementhr} -eq 0 ]; then
  echo "ERROR: increment hour (incrementhr) cannot be zero. Inifinite loop will result. Please modify. Exiting..."
  exit 1
fi

#----------------------------------------------------------------------
# End checks of user input
#----------------------------------------------------------------------

#----------------------------------------------------------------------
#  Begin work
#----------------------------------------------------------------------

# cd to working directory
cd ${DOMAINPATH}/postprd_${startdate}_F${fhr}_${startmin}
err1=$?
if test "$err1" -ne 0; then
    echo "ERROR: Could not 'cd' to working directory. Did you create directory: '${DOMAINPATH}/postprd'?  \
    Does '${DOMAINPATH}' exist?  Exiting... "
    exit 1
fi

# For GRIB2 the code reads a flat text tile to select variables for output.
#   The available fields are defined in post_avlbflds.xml -- while we
#   set a link to this file for reading during runtime it is not typical
#   for one to update this file, therefore the link goes back to the
#   program directory - this is true for params_grib2_tbl_new also - a
#   file which defines the GRIB2 table values
if [[ ${outFormat} == "grib2" ]]; then
   ln -fs ${txtCntrlFile} postxconfig-NT.txt
   ln -fs ${UPP_HOME}/parm/post_avblflds.xml post_avblflds.xml
   ln -fs ${UPP_HOME}/parm/params_grib2_tbl_new params_grib2_tbl_new
fi

# Link microphysics tables - code will use based on mp_physics option
# found in data
ln -fs ${UPP_HOME}/parm/nam_micro_lookup.dat .
ln -fs ${UPP_HOME}/parm/hires_micro_lookup.dat .

# link coefficients for crtm2 (simulated synthetic satellites)
CRTMDIR=${UPP_HOME}/crtm/fix
ln -fs $CRTMDIR/EmisCoeff/IR_Water/Big_Endian/Nalli.IRwater.EmisCoeff.bin           ./
ln -fs $CRTMDIR/EmisCoeff/MW_Water/Big_Endian/FASTEM4.MWwater.EmisCoeff.bin           ./
ln -fs $CRTMDIR/EmisCoeff/MW_Water/Big_Endian/FASTEM5.MWwater.EmisCoeff.bin           ./
ln -fs $CRTMDIR/EmisCoeff/MW_Water/Big_Endian/FASTEM6.MWwater.EmisCoeff.bin           ./
ln -fs $CRTMDIR/EmisCoeff/IR_Land/SEcategory/Big_Endian/NPOESS.IRland.EmisCoeff.bin           ./
ln -fs $CRTMDIR/EmisCoeff/IR_Snow/SEcategory/Big_Endian/NPOESS.IRsnow.EmisCoeff.bin           ./
ln -fs $CRTMDIR/EmisCoeff/IR_Ice/SEcategory/Big_Endian/NPOESS.IRice.EmisCoeff.bin           ./
ln -fs $CRTMDIR/AerosolCoeff/Big_Endian/AerosolCoeff.bin     ./
ln -fs $CRTMDIR/CloudCoeff/Big_Endian/CloudCoeff.bin         ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_g11.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g11.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_g12.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g12.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_g13.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g13.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_g15.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_g15.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_mt1r.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_mt1r.TauCoeff.bin    
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_mt2.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_mt2.TauCoeff.bin    
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/imgr_insat3d.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/imgr_insat3d.TauCoeff.bin    
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/amsre_aqua.SpcCoeff.bin  ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/amsre_aqua.TauCoeff.bin  ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/tmi_trmm.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/tmi_trmm.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmi_f13.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmi_f13.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmi_f14.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmi_f14.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmi_f15.SpcCoeff.bin    ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmi_f15.TauCoeff.bin    ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f16.SpcCoeff.bin   ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f16.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f17.SpcCoeff.bin   ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f17.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f18.SpcCoeff.bin   ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f18.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f19.SpcCoeff.bin   ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f19.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ssmis_f20.SpcCoeff.bin   ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ssmis_f20.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/seviri_m10.SpcCoeff.bin   ./   
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/seviri_m10.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/v.seviri_m10.SpcCoeff.bin   ./   
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/abi_gr.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/abi_gr.SpcCoeff.bin   ./
ln -fs $CRTMDIR/TauCoeff/ODPS/Big_Endian/ahi_himawari8.TauCoeff.bin   ./
ln -fs $CRTMDIR/SpcCoeff/Big_Endian/ahi_himawari8.SpcCoeff.bin   ./

#######################################################
# 1. Run UPP
#
# The UPP is used to read native GFS and LAM  model 
# output and put out isobaric state fields and derived fields.
#######################################################

export NEWDATE=$startdate

YYY=`echo $startdate | cut -c1-4`
MMM=`echo $startdate | cut -c5-6`
DDD=`echo $startdate | cut -c7-8`
HHH=`echo $startdate | cut -c9-10`

while [ $((10#${fhr})) -le $((10#${lastfhr})) ]; do

# Formatted fhr for filenames
fhr=`printf "%02i" ${fhr#0}`
fhour=`printf "%03i" ${fhr##0}`

NEWDATE=`date '+%Y%m%d%H' --date="$YYY$MMM$DDD $HHH $((10#${fhr})) hour"`

YY=`echo $NEWDATE | cut -c1-4`
MM=`echo $NEWDATE | cut -c5-6`
DD=`echo $NEWDATE | cut -c7-8`
HH=`echo $NEWDATE | cut -c9-10`
iHH=`echo $startdate | cut -c9-10`

echo 'NEWDATE' $NEWDATE
echo 'YY' $YY
export min=${startmin}

if [ `echo "${startdate}" | awk '/^[[:digit:]]{10}$/'` ]; then
    startdate=`echo "${startdate}" | sed 's/\([[:digit:]]\{2\}\)$/ \1/'`
  elif [ ! "`echo "${startdate}" | awk '/^[[:digit:]]{8}[[:blank:]]{1}[[:digit:]]{2}$/'`" ]; then
    ${ECHO} "ERROR: start time, '${startdate}', is not in 'yyyymmddhh' or 'yyyymmdd hh' format"
    exit 1
  fi
  START_TIME=`date -d "${startdate}"`

while [ $min -le $lastmin ]; do
timestr=`${DATE} +%Y-%m-%d_%H_%M_%S -d "${START_TIME}  ${fhr} hours ${min} minutes"`
timestr2=`${DATE} +%Y-%m-%d_%H:%M:%S -d "${START_TIME}  ${fhr} hours ${min} minutes"`

echo ${timestr}
# Create model file name (inFileName)
if [ ${model} == "GFS" ]; then
   if [ ${inFormat} == "netcdfpara" ]; then
      inFileName=${modelDataPath}/gfs.t00z.atmf${fhour}.nc
      flxFileName=${modelDataPath}/gfs.t00z.sfcf${fhour}.nc
   fi
elif [ ${model} == "LAM" ]; then
   if [ ${inFormat} == "netcdfpara" ]; then
      inFileName=${modelDataPath}/dynf${fhour}.nc
      flxFileName=${modelDataPath}/phyf${fhour}.nc
   fi
else
   inFileName=${modelDataPath}/wrfout_d01_${timestr}
   flxFileName=${modelDataPath}/wrfout_d01_${timestr}
fi

# Check if the files exist
if [[ ! -e ${inFileName} ]]; then
   echo "ERROR: Can't find 'inFileName': ${inFileName}. Directory or file does not exist.  Exiting..."
   echo "ERROR: Check if 'modelDataPath': ${modelDataPath} exists."
   exit 1
fi

if [[ ! -e ${flxFileName} ]]; then
   echo "ERROR: Can't find 'flxFileName': ${flxFileName}. Directory or file does not exist.  Exiting..."
   echo "ERROR: Check if 'modelDataPath': ${modelDataPath} exists."
   exit 1
fi

# Create itag based on user provided info. 
# Output format now set by user so if-block below uses this
# to generate the correct itag. 

if [[ ${outFormat} == "grib2" ]]; then
cat > itag <<EOF
&model_inputs
fileName='${inFileName}'
IOFORM='${inFormat}'
grib='${outFormat}'
DateStr='${YY}-${MM}-${DD}_${HH}:${min}:00'
MODELNAME='${tag}'
fileNameFlux='${flxFileName}'
fileNameFlat='postxconfig-NT.txt'
/
&nampgb
numx=${numx}
/
EOF
fi

#-----------------------------------------------------------------------
#   Run UPP.
#-----------------------------------------------------------------------

#----------------------------------------------------------------------
# There are two environment variables tmmark and COMSP
# RUN the upp.x executable. 
#----------------------------------------------------------------------

if [[ ${model} == "GFS" || ${model} == "LAM" || ${model} == "RAPR" ]]; then
     ${RUN_COMMAND} > upp.f${fhour}.out 2>&1
fi

# The prefixes are given in the postcntrl.xml file datset variable (GRIB2)

if [[ ${model} == "GFS" ]]; then
   mv GFSPRS.GrbF${fhr} GFSPRS.${fhour}
elif [ ${model} == "LAM" ]; then
   mv NATLEV${fhr}.tm00 NATLEV.${fhour}
   mv PRSLEV${fhr}.tm00 PRSLEV.${fhour}
elif [ ${model} == "RAPR" ]; then
   mv NATLEV${fhr}.tm00 NATLEV.${fhour}
   mv PRSLEV${fhr}.tm00 PRSLEV.${fhour}
fi

#
#----------------------------------------------------------------------
#   End of upp job
#----------------------------------------------------------------------

# check to make sure UPP was successful
if [[ ${model} == "GFS" ]]; then
    ls -l GFSPRS.${fhour}
    err1=$?
elif [ ${model} == "LAM" ]; then
    ls -l NATLEV.${fhour}
    err1=$? 
    ls -l PRSLEV.${fhour}
    err2=$? 
elif [ ${model} == "RAPR" ]; then
    ls -l WRFNAT*.tm00
    err1=$?
    ls -l WRFPRS*.tm00
    err2=$?    
fi

if [[ ${err1} -ne 0 || ${err2} -ne 0 ]]; then
    echo 'UPP FAILED, EXITTING'
    exit 1
else
    echo 'UPP SUCCEEDED'
# Add Subfix d01 to each GRIB file
    for file in ./WRF*; do
          if [ -f "$file" ]; then
              filename=$(basename -- "$file")
              new_filename="${filename}_d01"
              mv "$file" "$new_filename"
              echo "Renamed $filename to $new_filename"
          fi
      done

    echo $startdate
    initdate=$(echo "$startdate" | tr -d ' ')
    mkdir ${DOMAINPATH}/${initdate}
    mv WRF* ${DOMAINPATH}/${initdate}/
    rm -rf ${DOMAINPATH}/postprd_${initdate}_F${fhr}_${startmin}    
fi



let "min=min+$incrementmin"
done

fhr=$((10#${fhr}+$((${incrementhr}))))

NEWDATE=`date '+%Y%m%d%H' --date="$YYY$MMM$DDD $HHH $((10#${fhr})) hour"`

done

date
echo "End of Output Job" 
exit
