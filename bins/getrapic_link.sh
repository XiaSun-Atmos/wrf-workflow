#!/bin/bash
# Author: Xia Sun, xia.sun@noaa.gov
# Date: Oct, 2023
module load hpss
cd ${GFS_DIR}
echo ${GFS_DIR}
pwd
mkdir ${YYYY}${MM}${DD}${HH}00
cd ${YYYY}${MM}${DD}${HH}00
mkdir output

local_dir="/public/data/grids/rap/full/wrfnat/grib2"

year=${YYYY: -2}
echo $year
# Convert the input date to a Unix timestamp
input_date=$YYYY$MM$DD
timestamp=$(date -d "$input_date" +"%s")
echo $timestamp
# Convert the Unix timestamp to the Julian day of the year
julian_day=$(date -d "@$timestamp" +"%j")
echo ${julian_day}
# Derive prior Julian day of the year
timestamp2=$((timestamp-86400))
julian_day2=$(date -d "@$timestamp2" +"%j")

# Check if the file exists for 00 UTC with 21 hour forecast
for i in {00..21}; do
    file_path=${local_dir}/$year${julian_day}000000${i}
    echo ${file_path}
    wait_interval=10  # Adjust this to your preferred waiting interval in seconds
    echo "searching for ${i}."

    while [ ! -e "$file_path" ]; do
        echo "${i} fcst does not exist. Waiting for it to become available..."
        sleep $wait_interval
    done


    # Function to check if the file is available and not modified
    is_file_stable() {
        if [ -e "$file_path" ]; then
            initial_mtime=$(stat -c %Y "$file_path")  # Get the initial modification timestamp
            sleep 10  # Adjust the sleep duration to your desired interval
            current_mtime=$(stat -c %Y "$file_path")  # Get the current modification timestamp
            if [ "$initial_mtime" -eq "$current_mtime" ]; then
                return 0  # File is available and not modified
            fi
        fi
        return 1  # File is not available or has been modified
    }


    # Loop while the file is not stable
    while ! is_file_stable; do
        echo "File is not stable. Waiting for it to become stable..."
        sleep 10  # Adjust this to your desired interval
    done

    ln -s ${local_dir}/$year${julian_day}000000$i ${GFS_DIR}/${YYYY}${MM}${DD}${HH}00/output

done



# Check if the file exists for 21 UTC with 27-39 hour forecast
for i in {27..39}; do
    file_path=${local_dir}/$year${julian_day2}210000${i}
    echo ${file_path}
    wait_interval=10  # Adjust this to your preferred waiting interval in seconds
    echo "searching for ${i}."

    while [ ! -e "$file_path" ]; do
        echo "${i} fcst does not exist. Waiting for it to become available..."
        sleep $wait_interval
    done


    # Function to check if the file is available and not modified
    is_file_stable() {
        if [ -e "$file_path" ]; then
            initial_mtime=$(stat -c %Y "$file_path")  # Get the initial modification timestamp
            sleep 10  # Adjust the sleep duration to your desired interval
            current_mtime=$(stat -c %Y "$file_path")  # Get the current modification timestamp
            if [ "$initial_mtime" -eq "$current_mtime" ]; then
                return 0  # File is available and not modified
            fi
        fi
        return 1  # File is not available or has been modified
    }


    # Loop while the file is not stable
    while ! is_file_stable; do
        echo "File is not stable. Waiting for it to become stable..."
        sleep 10  # Adjust this to your desired interval
    done

    ln -s ${local_dir}/$year${julian_day2}210000$i ${GFS_DIR}/${YYYY}${MM}${DD}${HH}00/output

done




file_to_check="$year${julian_day}00000021"
# Check if the file exists
if [ -f "${GFS_DIR}/${YYYY}${MM}${DD}${HH}00/output/${file_to_check}" ]; then
    echo "File exists: $file_to_check"
    # Add your actions here if the file exists
else
    echo "Error: File does not exist: $file_to_check"
    exit 1  # Exit the script with an error code
fi
