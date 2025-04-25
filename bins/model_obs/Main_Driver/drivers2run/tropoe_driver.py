#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Driver for PSL TROPoe Information
-- Driver Details were Broken Up to Facilitate Data Analysis

###############################################################################
Created on Sun, Nov 12, 2023
@author: dmurray
"""

import os,sys
from datetime import date, timedelta
# find the config directory
sys.path.append(os.path.join(os.path.dirname(__file__),'../../config')) 
sys.path.append(os.path.join(os.path.dirname(__file__),'../../Utilities')) 
import config as cfg
###############################################################################
# Define the Predominant Working Directory 
CODEtree = getattr(cfg,'CODEdir')
DATAtree = getattr(cfg,'DATAdir')
MODELtree = getattr(cfg,'MODELdir')
TROPoetree = getattr(cfg,'TROPoedir')

###############################################################################
""" Denote Relevant Directories and Paths """
###############################################################################
sys.path.append(os.path.join(CODEtree,'config/parameter_files/'))
sys.path.append(os.path.join(CODEtree,'Main_Driver/general_functions/'))
sys.path.append(os.path.join(CODEtree,'Main_Driver/plotting_functions/'))
sys.path.append(os.path.join(CODEtree,'Utilities'))

for root,dirs,_ in os.walk(os.path.join(CODEtree,'Main_Driver/driver_instr_functions/')):
    if '__pycache__' in dirs: dirs.remove ('__pycache__') #do not add pycache directories
    for sdir in dirs:
        sys.path.append(os.path.join(root,sdir)) 

for root,dirs,_ in os.walk(os.path.join(CODEtree,'Main_Driver/plotting_functions/')):
    if '__pycache__' in dirs: dirs.remove ('__pycache__') #do not add pycache directories
    for sdir in dirs:
        sys.path.append(os.path.join(root,sdir)) 

for root,dirs,_ in os.walk(os.path.join(CODEtree,'Main_Driver/general_functions/')):
    if '__pycache__' in dirs: dirs.remove ('__pycache__') #do not add pycache directories
    for sdir in dirs:
        sys.path.append(os.path.join(root,sdir)) 

###############################################################################
# -> Image Output Directory (Server Based)
Sout = os.path.join(DATAtree,'Main_Driver/station_products/')
if not os.path.exists(Sout):
  os.makedirs(Sout)
#Smeas = os.path.join(DATAtree,'Measurement_Data/Stations/') # -> Datastream Output Directory
Smeas = TROPoetree # -> TROPoe data
Smodel = MODELtree # -> Model Input Directory
###############################################################################

###############################################################################
""" Load in Relevant Functions and Parameter Files """
import signal_param_list as signal_prm
import station_param_list as station_prm
######################################
from tropoe_main_driver import tropoe_driver

###############################################################################
""" Define Relevant Parameters for Main Driver Analyses """
###########################################################
# -> Initiate Analysis Period of Interest
PERoi = signal_prm.SGPper_oi; 
Tstrt = (date.today()-timedelta(days=PERoi)).strftime('%Y%m%d');
Tstop = (date.today()).strftime('%Y%m%d');
ANLYper = (Tstrt,Tstop);
###############################################################################
Doi = (Smeas,Smodel,Sout)

###############################################################################
""" Initiate TROPoe Data Analysis"""
for kind in ['ASSIST','MWR']:
   tropoe_driver(signal_prm.tropoe_oi,Doi,ANLYper,version=kind, overwrite=True)
###############################################################################
# Indicate that Analyses is Complete #
print('#####################################################################')
print(station_prm.EXPname+' TROPoe Data Analysis Complete: ' + Tstrt + '-' + Tstop)
print('#####################################################################')
