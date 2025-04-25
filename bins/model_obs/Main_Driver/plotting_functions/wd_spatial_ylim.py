#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Wind Direction Y-Limits Based on the Model 
and Observed Data

Do Not Need to Convert to U and V Wind Components Because of Correction to 0/
360 in the Colorbar Either Way

"""

import numpy as np
###############################################################################
RNDbse = 10;
TICKint = 5;

def wd_sptl_ylim(WDin):
    WDarr = [];
    ###########################################################################
    for w in WDin:
        WDarr = np.concatenate((WDarr,np.ravel(w)))
    ###########################################################################
    # Ensure Whether Several Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(WDarr)) > 1:
        # Determine the Maximum and Minimum Values of the Array
        #######################################################################
        MU = np.nanmean(WDarr);
        SIG = np.nanstd(WDarr);
        #######################################################################
        WDmax = np.ceil((MU + 3.1*SIG)/RNDbse)*RNDbse;
        WDmin = np.floor((MU - 3.1*SIG)/RNDbse)*RNDbse;        
        #######################################################################
        # Correct WDmin/WDmax Values
        WDmin = np.nanmax([0,WDmin])
        WDmax = np.nanmin([360,WDmax])
        #######################################################################
        #Determine Wind Direction Range
        WDrnge = WDmax - WDmin 
        #######################################################################
        stp = np.ceil(WDrnge/RNDbse)*(RNDbse/TICKint)
    else:
        WDmin = 0; WDmax = 360;
        WDrnge = WDmax - WDmin;
        #######################################################################
        stp = np.ceil(WDrnge/RNDbse)*(RNDbse/TICKint)

    WDtick = np.arange(WDmin,WDmax+stp,stp);    
    ###########################################################################
    return WDmin,WDmax,WDtick
