
#!/bin/bash


# for i in $(seq -w 00 48); do
#     ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFNAT${i}.tm00_d01 /data/input/wfip3.wrfnat.${i}.tm00_d01
#     ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFPRS${i}.tm00_d01 /data/input/wfip3.wrfprs.${i}.tm00_d01
#     ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFTWO${i}.tm00_d01 /data/input/wfip3.wrftwo.${i}.tm00_d01

#     python3 src/main.py -c ./conf -f /data/input/wfip3.wrfnat.${i}.tm00_d01 /data/input/wfip3.wrfprs.${i}.tm00_d01 /data/input/wfip3.wrftwo.${i}.tm00_d01
#     rm -rf /data/input/wfip3.wrfnat.${i}.tm00_d01 /data/input/wfip3.wrfprs.${i}.tm00_d01 /data/input/wfip3.wrftwo.${i}.tm00_d01
#     rm -rf /data/input/wfip3.wrfnat.${i}.tm00_d01*.idx /data/input/wfip3.wrfprs.${i}.tm00_d01*.idx /data/input/wfip3.wrftwo.${i}.tm00_d01*.idx 

# done

echo $SHELL


doy=$(date -d "${YYYY}-${MM}-${DD}" +%j)



    ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFNAThh.mm.tm00_d01 /data/input/wfip3km.wrfnat.hh.mm.tm00_d01
    ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFPRShh.mm.tm00_d01 /data/input/wfip3km.wrfprs.hh.mm.tm00_d01
    ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFTWOhh.mm.tm00_d01 /data/input/wfip3km.wrftwo.hh.mm.tm00_d01


    python3 src/main.py -c ./conf -f /data/input/wfip3km.wrfnat.hh.mm.tm00_d01 /data/input/wfip3km.wrfprs.hh.mm.tm00_d01 /data/input/wfip3km.wrftwo.hh.mm.tm00_d01
    rm -rf /data/input/wfip3km.wrfnat.hh.mm.tm00_d01 /data/input/wfip3km.wrfprs.hh.mm.tm00_d01 /data/input/wfip3km.wrftwo.hh.mm.tm00_d01
    rm -rf /data/input/wfip3km.wrfnat.hh.mm.tm00_d01*.idx /data/input/wfip3km.wrfprs.hh.mm.tm00_d01*.idx /data/input/wfip3km.wrftwo.hh.mm.15.tm00_d01*.idx
    # mv /data/output/hh_mm/profiler_sitesF.ncep_wfip3km.${YYYY}${doy}.00subhour.nc /data/output/profiler_sitesF.ncep_wfip3km.${YYYY}${doy}.${HH}hh.mm.nc
    # rm -rf /data/output/hh_mm

