#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Atmospheric Bias Y-Limits Based on the 
Model and Observed Data

"""

import numpy as np
###############################################################################

def atmo_blim(ATMOin,RNDbse,TICKint):
    ATMOarr = [];
    ###########################################################################
    for r in ATMOin:
        ATMOarr = np.concatenate((ATMOarr,np.ravel(r)))
    ###########################################################################
    # Ensure Whether Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(ATMOarr)) > 0:
        ATMOmax = np.ceil(np.nanmax(np.abs(ATMOarr))/RNDbse)*RNDbse;
        ATMOmin = -ATMOmax;
        #######################################################################
        """ # Determine the Maximum and Minimum Values of the Array
        ATMOmin = np.floor(np.nanmin(ATMOarr)/RNDbse)*RNDbse
        ATMOmax = np.ceil(np.nanmax(ATMOarr)/RNDbse)*RNDbse
        # No Corrections Needed for Min/Max Values (ie. Greater/Less than Some Value) """
        ###########################################################################
        #Determine the Relative Humidity Range
        ATMOrnge = ATMOmax - ATMOmin 
        #######################################################################
        stp = np.ceil(ATMOrnge/RNDbse)*(RNDbse/TICKint)
        #######################################################################
    else: #mainly for bias purposes
        ATMOmin = -4; ATMOmax = 4;
        ATMOrnge = ATMOmax - ATMOmin;
        #######################################################################
        stp = np.ceil(ATMOrnge/RNDbse)*(RNDbse/TICKint)

    ATMOtick = np.arange(ATMOmin,ATMOmax+stp,stp);    
    ###########################################################################
    return ATMOmin,ATMOmax,ATMOtick