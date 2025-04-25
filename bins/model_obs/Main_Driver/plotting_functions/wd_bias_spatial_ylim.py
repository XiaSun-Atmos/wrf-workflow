#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Wind Speed Bias Y-Limits Based on the 
Model and Observed Data
"""

import numpy as np
###############################################################################
RNDbse = 5;
TICKint = 5;

def wd_sptlbias_ylim(Bin):
    Barr = [];
    ###########################################################################
    for b in Bin:
        Barr = np.concatenate((Barr,np.ravel(b)))
    ###########################################################################
    # Ensure Whether Several Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(Barr)) > 1:
        # Determine the Maximum and Minimum Values of the Array
        #Bmin = np.floor(np.nanmin(Barr)/RNDbse)*RNDbse
        #Bmax = np.ceil(np.nanmax(Barr)/RNDbse)*RNDbse
        #######################################################################
        Bmin = np.floor((-2*np.nanstd(Barr))/RNDbse)*RNDbse;
        Bmax = -Bmin;
        #######################################################################
        #Determine Wind Speed Range
        Brnge = Bmax - Bmin 
        #######################################################################
        stp = np.ceil(Brnge/RNDbse)*(RNDbse/TICKint)
    else:
        Bmin = -15; Bmax = 15;
        Brnge = Bmax - Bmin;
        #######################################################################
        stp = np.ceil(Brnge/RNDbse)*(RNDbse/TICKint)

    Btick = np.arange(Bmin,Bmax+stp,stp)
    ###########################################################################
    return Bmin,Bmax,Btick