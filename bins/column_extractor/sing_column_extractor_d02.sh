
#!/bin/bash

    ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFNAThh.tm00_d02 /data/input/wfip1km.wrfnat.hh.tm00_d02
    ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFPRShh.tm00_d02 /data/input/wfip1km.wrfprs.hh.tm00_d02
    ln -s /upp_out/${YYYY}${MM}${DD}${HH}/WRFTWOhh.tm00_d02 /data/input/wfip1km.wrftwo.hh.tm00_d02

    python3 src/main.py -c ./conf -f /data/input/wfip1km.wrfnat.hh.tm00_d02 /data/input/wfip1km.wrfprs.hh.tm00_d02 /data/input/wfip1km.wrftwo.hh.tm00_d02
    rm -rf /data/input/wfip1km.wrfnat.hh.tm00_d02 /data/input/wfip1km.wrfprs.hh.tm00_d02 /data/input/wfip1km.wrftwo.hh.tm00_d02
    rm -rf /data/input/wfip1km.wrfnat.hh.tm00_d02*.idx /data/input/wfip1km.wrfprs.hh.tm00_d02*.idx /data/input/wfip1km.wrftwo.hh.tm00_d02*.idx 

    # ln -s /upp_out/2024080900/WRFNAT${i}.tm00_d02 /column_out/2024080900/d02/wfip3.wrfnat${i}.tm00_d02
    # ln -s /upp_out/2024080900/WRFPRS${i}.tm00_d02 /column_out/2024080900/d02/wfip3.wrfprs${i}.tm00_d02
    # ln -s /upp_out/2024080900/WRFTWO${i}.tm00_d02 /column_out/2024080900/d02/wfip3.wrftwo${i}.tm00_d02

    # python3 src/main.py -c ./conf -f /column_out/2024080900/d02/wfip3.wrfnat${i}.tm00_d02 /column_out/2024080900/d02/wfip3.wrfprs${i}.tm00_d02 /column_out/2024080900/d02/wfip3.wrftwo${i}.tm00_d02 -l /column_out/2024080900/d02 -o /column_out/2024080900/d02

