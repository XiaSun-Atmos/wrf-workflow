#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract Relevant PSL HRRR Model Data
###############################################################################
Created on Fri Dec 13 14:07:44 2019
@author: jduncan
"""
import tarfile, os
import numpy as np
from datetime import date
#####################################
from pslmdl_main_nc import pslmdl_nc

def pslmdl_ingest(dir_load,avail_flist,doi):
    #input includes path to mdl files (dir_load), available file list (avail_flist),
    #and the date of interedst (doi)    
    ###########################################################################
    # Convert Input Date of Interest to Day of Year
    dt_oi = date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8]))
    doy = dt_oi.strftime('%j')
    year = doi[:4]
    ########################################################
    # Define File of Interest
    #foi = doi[0:4]+doy
    foi = doi
    fil2proc = [fname for fname in avail_flist if foi in fname]
    ###########################################################################
    if len(fil2proc) == 0:
        """ Pre-Define Empty (i.e. length = 0) Tuples Herein """
        ########################################################
        PSLini = [];
        PSLxy = (); PSLloc = ();
        PSLpr = ();

    else:
        # Sort to Ensure Data Files Correctly Concatenated #
        #### ftime = [hfil.split('.')[3] for hfil in Thr_nc]
        #### Thr_nc = [x for _,x in sorted(zip(np.asarray(ftime),Thr_nc))]
        ftime = [hfil.split('.')[3] for hfil in fil2proc]
        Thr_nc = [x for _,x in sorted(zip(np.asarray(ftime),fil2proc))]
        
        """ Extract Data from Embedded .nc Files """
        #######################################################################
        # Pre-establish Tuples of Interest #
        PSLini = [];
        PSLxy = (); PSLloc = ();
        PSLpr = ();
        for hnc in Thr_nc:
            # Extract netCDF to the Local Working Directory #
            #Toi.extract(hnc,Ttemp)
            # Process .nc File
            NCxy, NCloc, NCpr, = pslmdl_nc(os.path.join(dir_load,year),hnc)
            # Remove .nc File
            #os.remove(os.path.join(Ttemp,hnc))
            ###################################################################
            """ Create Large Scale Tuple for Each Data File """
            ###################################################################
            #-> Extract Initialization Time 
            PSLini.append(hnc.split('.')[-2][0:2])
            ###################################################################
            PSLxy = PSLxy + (NCxy,)
            PSLloc = PSLloc + (NCloc,)
            PSLpr = PSLpr + (NCpr,)
            #commas are needed
        #######################################################################
        #Toi.close();
        
    """ Return Tuples of Interest """
    return PSLini,PSLxy, PSLloc, PSLpr

