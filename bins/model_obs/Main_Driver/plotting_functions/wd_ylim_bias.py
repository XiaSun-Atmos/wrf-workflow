#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Wind Direction Bias Y-Limits Based on the 
Model and Observed Data
"""

import numpy as np
###############################################################################
RNDbse = 5;
TICKint = 5;

def wd_blim(Bin):
    Barr = [];
    ###########################################################################
    for b in Bin:
        Barr = np.concatenate((Barr,np.ravel(b)))
    ###########################################################################
    # Ensure Whether Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(Barr)) > 0:
        Bmax = np.ceil(np.nanmax(np.abs(Barr))/RNDbse)*RNDbse;
        Bmin = -Bmax;
        #######################################################################
        """# Determine the Maximum and Minimum Values of the Array
        Bmin = np.floor(np.nanmin(Barr)/RNDbse)*RNDbse
        Bmax = np.ceil(np.nanmax(Barr)/RNDbse)*RNDbse """
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

    # Check for the case where there is only one value
    if (stp == 0):
        Bmin = -15; Bmax = 15;
        Brnge = Bmax - Bmin;
        stp = np.ceil(Brnge/RNDbse)*(RNDbse/TICKint)

    Btick = np.arange(Bmin,Bmax+stp,stp)
    ###########################################################################
    return Bmin,Bmax,Btick
