#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract Relevant ESRL Model Data
###############################################################################
Created on Fri Dec 13 14:07:44 2019
@author: jduncan
"""
import tarfile, os
import numpy as np
from datetime import date
#####################################
from esrl_main_nc import esrl_nc

def esrl_ingest(dir_load,avail_flist,doi):
    #input includes path to mdl files (dir_load), available file list (avail_flist),
    #and the date of interedst (doi)    
    ###########################################################################
    # Convert Input Date of Interest to Day of Year
    dt_oi = date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8]))
    doy = dt_oi.strftime('%j')
    year = doi[:4]
    ########################################################
    # Define File of Interest
    foi = doi[0:4]+doy
    fil2proc = [fname for fname in avail_flist if ((foi in fname) & (fname.endswith(".nc")))]
    ###########################################################################
    if len(fil2proc) == 0:
        """ Pre-Define Empty (i.e. length = 0) Tuples Herein """
        ########################################################
        ESRLini = [];
        ESRLxy = (); ESRLloc = ();
        ESRLtpmx = (); ESRLwind = (); ESRLflux = (); ESRLrad = ();

    else:
        """ Establish Temporary Working Directory to Unpack and Analyze """
        #### Ttemp = os.path.join(dir_load,'temp'); 
        #### if not os.path.exists(Ttemp):
        ####     os.makedirs(Ttemp)
        #######################################################################
        #### Toi = tarfile.open(os.path.join(dir_load,fil2proc[0]),'r:gz');
        # -- extract the names of the embedded .nc files -- #
        #### Thr_nc = Toi.getnames();
        #######################################################################
        # Sort to Ensure Data Files Correctly Concatenated #
        #### ftime = [hfil.split('.')[3] for hfil in Thr_nc]
        #### Thr_nc = [x for _,x in sorted(zip(np.asarray(ftime),Thr_nc))]
        ftime = [hfil.split('.')[3] for hfil in fil2proc]
        Thr_nc = [x for _,x in sorted(zip(np.asarray(ftime),fil2proc))]
        
        """ Extract Data from Embedded .nc Files """
        #######################################################################
        # Pre-establish Tuples of Interest #
        ESRLini = [];
        ESRLxy = (); ESRLloc = ();
        ESRLtpmx = (); ESRLwind = (); ESRLflux = (); ESRLrad = ();
        for hnc in Thr_nc:
            # Extract netCDF to the Local Working Directory #
            #Toi.extract(hnc,Ttemp)
            # Process .nc File
            #print('processing',hnc)
            NCxy, NCloc, NCtpmx, NCwind, NCflux, NCrad = esrl_nc(os.path.join(dir_load,year),hnc)
            # Remove .nc File
            #os.remove(os.path.join(Ttemp,hnc))
            ###################################################################
            """ Create Large Scale Tuple for Each Data File """
            ###################################################################
            #-> Extract Initialization Time 
            ESRLini.append(hnc.split('.')[-2][0:2])
            ###################################################################
            ESRLxy = ESRLxy + (NCxy,)
            ESRLloc = ESRLloc + (NCloc,)
            ESRLtpmx = ESRLtpmx + (NCtpmx,)
            ESRLwind = ESRLwind + (NCwind,)
            ESRLflux = ESRLflux + (NCflux,)
            ESRLrad = ESRLrad +(NCrad,)
            #commas are needed
        #######################################################################
        #Toi.close();
        
    """ Return Tuples of Interest """
    return ESRLini,ESRLxy, ESRLloc, ESRLtpmx, ESRLwind, ESRLflux, ESRLrad

