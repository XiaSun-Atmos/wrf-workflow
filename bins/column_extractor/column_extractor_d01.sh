#!/bin/bash

export column_out="/mnt/lfs5/BMC/rtwbl/Xia.Sun/wfip3_cases/outputs/column_out"
export upp_out="/mnt/lfs5/BMC/rtwbl/Xia.Sun/wfip3_cases/outputs/upp_out"
export column_extract="/mnt/lfs5/BMC/rtwbl/Xia.Sun/wfip3_cases//workflow/bins/column_extractor"
export SINGULARITYENV_APPEND_PATH="/opt/conda/bin:${column_extract}"
# export SINGULARITY_BIND="${upp_out}:/upp_out,${column_out}/${YYYY}${MM}${DD}${HH}/d01:/data"
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/input
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/output

module load cdo
module load nco

# cd ${column_extract}
# singularity exec -H ${column_extract} ${column_extract}/column_extractor_subh_latest.sif ${column_extract}/sing_column_extractor_d01.sh

column_extract_subh (){

cd ${column_extract}
cp -r ${column_extract}/sing_column_extractor_subh_d01.sh ${column_extract}/sing_column_extractor_d01_${MM}${DD}${1}_${2}.sh

sed -i "s/hh/${1}/g" sing_column_extractor_d01_${MM}${DD}${1}_${2}.sh
sed -i "s/mm/${2}/g" sing_column_extractor_d01_${MM}${DD}${1}_${2}.sh

mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_${2}
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_${2}/input
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_${2}/output

export SINGULARITY_BIND="${upp_out}:/upp_out,${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_${2}:/data"

singularity exec -H ${column_extract} ${column_extract}/column_extractor_subh_latest.sif ${column_extract}/sing_column_extractor_d01_${MM}${DD}${1}_${2}.sh

subh=$(echo "scale=2; ${1} + ${2} / 60" | bc)

cd ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_${2}/output

column_file=$(ls *.nc)

echo ${column_file}

ncap2 -h -s 'forecast=float(forecast)' ${column_file} -O ${column_file}

ncap2 -h -s "forecast(:)=${subh}" ${column_file} -O ${column_file}

new_file="${column_file/subhour/${1}.${2}}"
mv ${column_file} ${column_out}/${YYYY}${MM}${DD}${HH}/d01/output/${new_file}

rm -rf ${column_extract}/sing_column_extractor_d01_${MM}${DD}${1}_${2}.sh
rm -rf  ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_${2}



}

export -f column_extract_subh

parallel -u column_extract_subh ::: {00..47} ::: 15 30 45


column_extract (){

cd ${column_extract}
cp -r ${column_extract}/sing_column_extractor_d01.sh ${column_extract}/sing_column_extractor_d01_${MM}${DD}${1}_00.sh

sed -i "s/hh/${1}/g" sing_column_extractor_d01_${MM}${DD}${1}_00.sh

mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_00
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_00/input
mkdir ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_00/output

export SINGULARITY_BIND="${upp_out}:/upp_out,${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_00:/data"

singularity exec -H ${column_extract} ${column_extract}/column_extractor_subh_latest.sif ${column_extract}/sing_column_extractor_d01_${MM}${DD}${1}_00.sh

subh=$(echo "scale=2; ${1} + 0 / 60" | bc)

cd ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_00/output

column_file=$(ls *.nc)

echo ${column_file}

ncap2 -h -s 'forecast=float(forecast)' ${column_file} -O ${column_file}

ncap2 -h -s "forecast(:)=${subh}" ${column_file} -O ${column_file}

new_file="${column_file/subhour/${1}.00}"
mv ${column_file} ${column_out}/${YYYY}${MM}${DD}${HH}/d01/output/${new_file}

rm -rf ${column_extract}/sing_column_extractor_d01_${MM}${DD}${1}_00.sh

rm -rf  ${column_out}/${YYYY}${MM}${DD}${HH}/d01/${1}_00

}

export -f column_extract

parallel -u column_extract ::: {00..48}