#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Purpose of this Function is to Ingest Spatial Meaurement Data (e.g. 
Profiling Lidar Data) and Convert it to Hourly Data to Enable Comparison between
 the Measured and Model Data
###############################################################################
-> This Function Will Focus on Wind Direction Data (2d) as it Cannot be Readily 
-> Averaged
###############################################################################
Created on Wed Dec 18 14:47:46 2019
@author: jduncan

This Script Sets up the Inequality to Handle Those Measurements whose Time 
Stamp Represents the 'End' of the Averaging Interval (i.e. Radar Wind Profiler)
or the Doppler Lidar who Outputs at 15-min and whose Time Stamp Represents the 
Middle of the Averaging Period. The Lidar Time stamp never falls exactly on a 
one-hour interval and therefore the inequality employed below is appropriate for
both profiling datastreams.

The Employed Inequality is:
    FCSTtime - 1800s < MEAStime <= FCSTtime + 1800s

"""

import numpy as np
import warnings
############################
from met_to_uv import wd2uv
from uv_to_met import uv2wd
############################
HRsec = 60*60;
AREAreq = 75; 

def wd_meas2grid(MEASwd_var,MEASgrd,MDLgrd):
    """ Define Function Input
    ###########################################################################
    MEASwd_var = Gridded Spatial Data (WD) 
    MEASgrd = Tuple Depicting Measurement Grid (Hgt,Tme) Arrays
    MDLgrd = Tuple Decpiting Model Grid (Hgt,Tme) Arrays
    ########################################################################"""
    # Unpack Input Tuples
    MDLhgt = MDLgrd[0]; MDLtme = MDLgrd[1]; # modeled coordinates
    #############################################################
    # reference measured meshgrid
    MEASh_mesh = MEASgrd[0]; MEASt_mesh = MEASgrd[1]; # measured coordinates
    ###########################################################################
    # Preallocate Averaging Array 
    MEASwd_mu = np.ones((len(MDLhgt),len(MDLtme)))*np.nan

    ########################################################################"""
    """ Decompose Wind Direction Data into U and V Wind Components """
    WSscale = np.ones(MEASwd_var.shape); #array of ones (required for scaler averaging)
    Uvec,Vvec = wd2uv(WSscale,MEASwd_var)
    
    ###########################################################################
    #-> preallocate averaging array
    MEASwd_mu = np.ones((len(MDLhgt),len(MDLtme)))*np.nan
    
    ###########################################################################
    for tod in range(0,len(MDLtme)): #for each time of day
        """ Define Relevant Time Averaging Period """
        TMEind_ref = np.logical_and(MEASt_mesh > MDLtme[tod] - HRsec/2,MEASt_mesh <= MDLtme[tod] + HRsec/2)
        
        """ Define Relevant Measurement Heights for Averaging """
        for hog in range(0,len(MDLhgt)): #for each height of grid
            ###################################################################
            # Define Half-Way Points Between Neighboring Bins
            if hog != 0 and hog != len(MDLhgt)-1:
                LWhgt_bin = (MDLhgt[hog] - MDLhgt[hog-1])/2;
                UPhgt_bin = (MDLhgt[hog+1] - MDLhgt[hog])/2;
                ###############################################################
                # Determine Relevant Averaging Reference Indices
                HGTind_ref = np.logical_and(MEASh_mesh >= MDLhgt[hog] - LWhgt_bin,\
                                            MEASh_mesh < MDLhgt[hog] + UPhgt_bin)
            else:
                if hog == 0:
                    # Simply Use the Upper Height Delta
                    HGTbin = (MDLhgt[hog+1] - MDLhgt[hog])/2; 
                elif hog == len(MDLhgt)-1:
                    # Simply Use the Lower Height Delta
                    HGTbin = (MDLhgt[hog] - MDLhgt[hog-1])/2;     
                ###############################################################
                # Determine Relevant Averaging Reference Indices
                HGTind_ref = np.logical_and(MEASh_mesh >= MDLhgt[hog] - HGTbin,\
                                            MEASh_mesh < MDLhgt[hog] + HGTbin)
            ###################################################################
            """ Define the Appropriate Atmospheric Variable Average """
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',category=RuntimeWarning)
                ###########################################################
                # create reference array
                MEASin_grd = MEASwd_var[np.logical_and(TMEind_ref,HGTind_ref)];
                ###############################################################
                if len(MEASin_grd)>0: #require at least one observation
                    # determine number of defined indices
                    MEASin_def = np.count_nonzero(~np.isnan(MEASin_grd))
                    MEASin_pct = MEASin_def/len(MEASin_grd)*100
                    # require a certain percent of non-nan measurements in area
                    if MEASin_pct >= AREAreq:
                        Umu = np.nanmean(Uvec[np.logical_and(TMEind_ref,HGTind_ref)]);
                        Vmu = np.nanmean(Vvec[np.logical_and(TMEind_ref,HGTind_ref)]);
                        #######################################################
                        MEASwd_mu[hog,tod] = uv2wd(Umu,Vmu);
                        
    """ Return the Modified Grid """
    return MEASwd_mu
