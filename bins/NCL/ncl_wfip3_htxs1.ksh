#!/bin/ksh --login

# Load modules
export MODULE_FILE="/home/role.rtrr/PARM_EXEC/modulefiles/modulefile.jet.NCL"
source ${MODULE_FILE}

# Make sure we are using GMT time zone for time computations
export TZ="GMT"
# export NCARG_ROOT="/apps/ncl/6.5.0-CentOS6.10_64bit_nodap_gnu447"
# export NCARG_LIB="/apps/ncl/6.5.0-CentOS6.10_64bit_nodap_gnu447/lib"
export NCL_HOME="/whome/Brian.D.Jamison/fim/svncode/ncl/fimall"
#export NCL_HOME="/lfs5/BMC/rtwbl/WFIP3/workflow/bins/fimall"
# export UDUNITS2_XML_PATH=$NCARG_ROOT/lib/ncarg/udunits/udunits2.xml

# Set up paths to shell commands
LS=/bin/ls
LN=/bin/ln
RM=/bin/rm
MKDIR=/bin/mkdir
CP=/bin/cp
MV=/bin/mv
ECHO=/bin/echo
CAT=/bin/cat
GREP=/bin/grep
CUT=/bin/cut
AWK="/bin/gawk --posix"
SED=/bin/sed
DATE=/bin/date
BC=/usr/bin/bc
NCL=`which ncl`
CTRANS=`which ctrans`
PS2PDF=/usr/bin/ps2pdf
CONVERT=`which convert`
MONTAGE=`which montage`
# PATH=${NCARG_ROOT}/bin:${PATH}

# typeset -RZ2 FCST_TIME

# ulimit -s 512000
ulimit -s 1024000

# Settings for testing
# NATINFILEPATH=/lfs1/BMC/wrfruc/jamison/adb_graphics/testwfip
# WORKDIR_ROOT=/lfs1/BMC/nrtrr/ncl_graphics/testwfipworkdir
# OUTPUT_ROOT=/lfs1/BMC/nrtrr/ncl_graphics/wfipxs
# EXE_ROOT=/misc/whome/role.wrfruc/bin/ncl/nclwfip
# START_TIME=2024033000
# FCST_TIME=00

# Print run parameters
${ECHO}
${ECHO} "ncl.ksh started at `${DATE}`"
${ECHO}

# Check to make sure the EXE_ROOT var was specified
if [ ! -d ${EXE_ROOT} ]; then
  ${ECHO} "ERROR: EXE_ROOT, '${EXE_ROOT}', does not exist"
  exit 1
fi

# Check to make sure that the DATAHOME exists
if [ ! -d ${NATINFILEPATH} ]; then
  ${ECHO} "ERROR: NATINFILEPATH, '${NATINFILEPATH}', does not exist"
  exit 1
fi
# If START_TIME is not defined, use the current time
if [ ! "${START_TIME}" ]; then
  ${ECHO} "START_TIME not defined - get from date"
  START_TIME=$( date +"%Y%m%d %H" )
  START_TIME=$( date +"%Y%m%d%H" -d "${START_TIME}" )
else
  ${ECHO} "START_TIME defined and is ${START_TIME}"
  START_TIME=$( date +"%Y%m%d %H" -d "${START_TIME%??} ${START_TIME#????????}" )
  START_TIME=$( date +"%Y%m%d%H" -d "${START_TIME}" )
fi

# Print out times
# ${ECHO} "   START TIME = "`${DATE} +%Y%m%d%H -d "${START_TIME}"`
${ECHO} "   START_TIME = ${START_TIME}"
${ECHO} "   FCST_TIME = ${FCST_TIME}"

# Set up the work directory and cd into it
workdir=${WORKDIR_ROOT}/${START_TIME}${FCST_TIME}htxs1
if [ -d ${workdir} ]; then
  ${RM} -rf ${workdir}
fi
${MKDIR} -p ${workdir}
cd ${workdir}

curdir=`pwd`
if [ "${curdir}" != "${workdir}" ]; then
  ${ECHO} "curdir = ${curdir}"
  ${ECHO} "workdir = ${workdir}"
  exit
fi

# Link to input file
FILE_NAME=WRFNAT${FCST_TIME}.tm00_d01
if [ ! -f ${NATINFILEPATH}/${START_TIME}/${FILE_NAME} ]; then
  ${ECHO} "ERROR: '${NATINFILEPATH}/${START_TIME}/${FILE_NAME}', does not exist"
  ${RM} -rf ${workdir}
  exit
fi  
if [ -f nat.grb ]; then
  ${RM} nat.grb
fi  

${LN} -s ${NATINFILEPATH}/${START_TIME}/${FILE_NAME} nat.grb
${ECHO} "${workdir}/nat.grb" > nat_file.txt

${ECHO} "NATINFILEPATH = ${NATINFILEPATH}"
${ECHO} "FILE_NAME = ${FILE_NAME}"

${LN} -s ${EXE_ROOT}/WRFUserARW.ncl WRFUserARW.ncl
${LN} -s ${EXE_ROOT}/stdatm.txt stdatm.txt

set -A ncgms sitesA_htxswind  \
             sitesA_htxscond  \
             sitesA_htxsclwmr \
             sitesA_htxsrwmr  \
             sitesA_htxssnmr  \
             sitesA_htxsgrmr  \
             sitesA_htxssmoke

set -A pngs sitesA_htxswind.000001.png \
            sitesA_htxswind.000002.png \
            sitesA_htxswind.000003.png \
            sitesA_htxscond.000001.png \
            sitesA_htxscond.000002.png \
            sitesA_htxscond.000003.png \
            sitesA_htxsclwmr.000001.png \
            sitesA_htxsclwmr.000002.png \
            sitesA_htxsclwmr.000003.png \
            sitesA_htxsrwmr.000001.png \
            sitesA_htxsrwmr.000002.png \
            sitesA_htxsrwmr.000003.png \
            sitesA_htxssnmr.000001.png \
            sitesA_htxssnmr.000002.png \
            sitesA_htxssnmr.000003.png \
            sitesA_htxsgrmr.000001.png \
            sitesA_htxsgrmr.000002.png \
            sitesA_htxsgrmr.000003.png \
            sitesA_htxssmoke.000001.png \
            sitesA_htxssmoke.000002.png \
            sitesA_htxssmoke.000003.png

set -A webnames htxswind_wf1A  \
                htxswind_wf2A  \
                htxswind_wf3A  \
                htxscond_wf1A  \
                htxscond_wf2A  \
                htxscond_wf3A  \
                htxsclwmr_wf1A  \
                htxsclwmr_wf2A  \
                htxsclwmr_wf3A  \
                htxsrwmr_wf1A  \
                htxsrwmr_wf2A  \
                htxsrwmr_wf3A  \
                htxssnmr_wf1A  \
                htxssnmr_wf2A  \
                htxssnmr_wf3A  \
                htxsgrmr_wf1A  \
                htxsgrmr_wf2A  \
                htxsgrmr_wf3A  \
                htxssmoke_wf1A  \
                htxssmoke_wf2A  \
                htxssmoke_wf3A

ncl_error=0

# Run the NCL scripts for each plot
cp ${EXE_ROOT}/names_grib2.txt .
cp ${EXE_ROOT}/parms*.txt .

i=0
while [ ${i} -lt ${#ncgms[@]} ]; do

  plot=${ncgms[${i}]}
  ${ECHO} "Starting rr_${plot}.ncl at `${DATE}`"
  ${NCL} < ${EXE_ROOT}/rr_${plot}.ncl
  error=$?
  if [ ${error} -ne 0 ]; then
    ${ECHO} "ERROR: rr_${plot} crashed!  Exit status=${error}"
    ncl_error=${error}
  fi
  ${ECHO} "Finished rr_${plot}.ncl at `${DATE}`"

  (( i=i + 1 ))

done

# Copy png files to their proper names
if [ ! -d ${OUTPUT_ROOT}/htxs ]; then
  ${MKDIR} -p ${OUTPUT_ROOT}/htxs
  ${CHMOD} 755 ${OUTPUT_ROOT}/htxs
fi  

i=0 
while [ ${i} -lt ${#pngs[@]} ]; do
  pngfile=${pngs[${i}]}
  webfile=${OUTPUT_ROOT}/htxs/${webnames[${i}]}_f0${FCST_TIME}.png
#  webfile=${webnames[${i}]}_f${FCST_TIME}.png    # for testing
  ${MV} ${pngfile} ${webfile}

  (( i=i + 1 ))
done

# Remove the workdir
${RM} -rf ${workdir}

${ECHO} "ncl.ksh completed at `${DATE}`"

exit ${ncl_error}
