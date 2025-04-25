#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Date: June, 2024
# Purpurse: retrieve RAP IC from HPSS
module load hpss

date=${YYYY}${MM}${DD}${HH}
forecast=${fcsthrs}

# Extract the year (first four characters)
year=${date:0:4}
mm=${date:4:2}
dd=${date:6:2}
yy=${date:2:2}
hh=${date:8:2}
echo $year
echo $mm
echo $dd
echo $yy
echo $hh

dir=${GFS_DIR}
mkdir ${dir}/${date}

cd ${dir}/${date}

export date
export year
export mm
export dd
export hh


hsi_get (){

    enddate=$(date -d "${year}-${mm}-${dd} ${hh}:00:00 ${1} hours" +"%Y%m%d%H")
    echo ${enddate}
    zipfile=${enddate:0:4}${enddate:4:2}${enddate:6:2}${enddate:8:2}00.zip

    timestamp=$(date -d "${enddate:0:8}" +"%s")
    julian_day=$(date -d "@$timestamp" +"%j")
    endhh=${enddate:8:2}
    echo ${julian_day}
    filename="${enddate:2:2}${julian_day}${endhh}000000"
    echo ${filename}

    if [ -f "$filename" ]; then
    echo "$filename exists, skipping..."
    else

    hsi get /BMC/fdr/Permanent/${enddate:0:4}/${enddate:4:2}/${enddate:6:2}/grib/ftp_rap_hyb/7/0/105/0_794802_32769/${zipfile}


    unzip -j "$zipfile" "$filename" -d "./"
    rm -rf $zipfile
    fi
}

export -f hsi_get
parallel -u hsi_get ::: {00..48}




# hsi get /BMC/fdr/Permanent/${year}/${mm}/${dd}/grib/ftp_rap_hyb/7/0/105/0_794802_32769/${year}${mm}${dd}0000.zip
# hsi get /BMC/fdr/Permanent/${year}/${mm}/${dd}/grib/ftp_rap_hyb/7/0/105/0_794802_32769/${year}${mm}${dd}0300.zip
# # Check the exit status of the hsi get command
# if [ $? -eq 0 ]; then
#     echo "Success: File retrieval was successful."
# else
#     echo "Error: File retrieval failed. Exit status: $?"
#     exit 1  # Exit the script with an error code
# fi

# zipfile=${year}${mm}${dd}0000.zip
# for i in $(seq -w 0 21); do
#           filename="$yy${julian_day}000000$i"
#           echo $filename
#       	rm -r $yy${julian_day}000000$i
#           unzip -j "$zipfile" "$filename" -d "./"
#           # Check the exit status of gunzip
#                 if [ $? -eq 0 ]; then
#                     echo "Success: File is uncompressed correctly."
#                 else
#                         rm -r $yy${julian_day}000000$i
#                         unzip -j "$zipfile" "$filename" -d "./"
#                     echo "redo the uncompressing"
#                     if [ $? -eq 1 ]; then
#                         echo "Failed: File is uncompressed incorrectly."
#                         fi

#                 fi
#           echo $i
# done


# zipfile2=${year}${mm}${dd}0300.zip
# for i in $(seq -w 21 31); do
#           filename="$yy${julian_day}030000$i"
#           echo $filename
#       	rm -r $yy${julian_day}030000$i
#           unzip -j "$zipfile2" "$filename" -d "./"
#           # Check the exit status of gunzip
#                 if [ $? -eq 0 ]; then
#                     echo "Success: File is uncompressed correctly."
#                 else
#                         rm -r $yy${julian_day}030000$i
#                         unzip -j "$zipfile" "$filename" -d "./"
#                     echo "redo the uncompressing"
#                     if [ $? -eq 1 ]; then
#                         echo "Failed: File is uncompressed incorrectly."
#                         fi

#                 fi
#           echo $i
# done

# rm -rf ${year}${mm}${dd}0000.zip
# rm -rf ${year}${mm}${dd}0300.zip
