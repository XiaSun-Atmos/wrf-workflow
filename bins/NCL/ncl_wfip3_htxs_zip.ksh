#!/bin/ksh --login

DATE=/bin/date
ECHO=/bin/echo

ulimit -s 512000

# Print run parameters
${ECHO}
${ECHO} "ncl_wfip3_htxs_zip.ksh started at `${DATE}`"
${ECHO}
${ECHO} "OUTPUT_ROOT = ${OUTPUT_ROOT}"

# Check to make sure that the OUTPUT_ROOT exists
if [ ! -d ${OUTPUT_ROOT} ]; then
  ${ECHO} "ERROR: OUTPUT_ROOT, '${OUTPUT_ROOT}', does not exist"
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
${ECHO} "   START_TIME = ${START_TIME}"

set -A domains htxs

zip_error=0

# zip up the files in each domain

i=0
while [ ${i} -lt ${#domains[@]} ]; do
  dir=${OUTPUT_ROOT}/${START_TIME}/nclprd_d01/${domains[${i}]}
  if [ -d ${dir} ]; then
    cd ${dir}
    if (( `ls *.png 2> /dev/null|wc -l` ));then
      zip -n .png files.zip * -i \*.png 
      zip_error=$?
      if [ zip_error -ne 0 ]; then
        ${ECHO} "ERROR - zip failed!"
        ${ECHO} " zip_error = ${zip_error}"
      else
        ${ECHO} "SUCCESS - zip file created"
        ${ECHO} " zip_error = ${zip_error}"
        rm -f *.png
      fi
    else
      ${ECHO} "no files to zip -- exiting"
    fi 
  else
    ${ECHO} "${dir} does not exist"
  fi 

  (( i=i + 1 ))
done

${ECHO} "ncl_wfip3_htxs_zip.ksh completed at `${DATE}`"

exit ${zip_error}
