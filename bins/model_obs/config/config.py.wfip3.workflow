"""
# Defining Relevant Top Level Directories
###############################################################################
Created on Wed Feb 14 2024
@author: dmurray

Adapted for usange in workflow Fri Dec 6 2024
@author: xsun
"""
##
##  Set the paths to the appropriate data/files, then copy to config.py
##

import os

bin_dir = os.environ.get('BIN_DIR')
obs_dir = os.environ.get('OBS_DIR')
out_dir = os.environ.get('OUT_DIR')

CODEdir = '{bin_dir}/model_obs'
DATAdir = '{obs_dir}'
MODELdir = '{out_dir}/column_out'

TROPoedir = '/Projects/WFIP3/processed/realtime/TROPoe/'
