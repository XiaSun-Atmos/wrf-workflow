#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Wind Direction Y-Limits Based on the Model 
and Observed Data

Y-Lims Will Be Based on the Assumption that there are Two Subplots 
--> This Can Potentially be Optimized to Allow More Y-ticks
"""

import numpy as np
###############################################################################
RNDbse = 10;
TICKint = 5;

def wd_ylim(WDin):
    WDarr = [];
    ###########################################################################
    for w in WDin:
        WDarr = np.concatenate((WDarr,np.ravel(w)))
    ###########################################################################
    # Ensure Whether Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(WDarr)) > 0:
        # Determine the Maximum and Minimum Values of the Array
        WDmin = np.floor(np.nanmin(WDarr)/RNDbse)*RNDbse
        WDmax = np.ceil(np.nanmax(WDarr)/RNDbse)*RNDbse
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

    # Check for the case where there is only one value
    if (stp == 0):
        WDmin = 0; WDmax = 360;
        WDrnge = WDmax - WDmin;
        #######################################################################
        stp = np.ceil(WDrnge/RNDbse)*(RNDbse/TICKint)

    WDtick = np.arange(WDmin,WDmax+stp,stp);    
    ###########################################################################
    return WDmin,WDmax,WDtick
