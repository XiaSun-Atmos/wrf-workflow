#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Limits Based on the Observed Atmospheric
Bias Data
"""

import numpy as np
###############################################################################

def atmo_sptlbias_ylim(ATMOin,RNDbse,TICKint):
    ATMOarr = [];
    ###########################################################################
    for r in ATMOin:
        ATMOarr = np.concatenate((ATMOarr,np.ravel(r)))
    ###########################################################################
    # Ensure Whether Several Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(ATMOarr)) > 1:
        # Determine the Maximum and Minimum Values of the Array
        ATMOmin = np.floor((-3*np.nanstd(ATMOarr))/RNDbse)*RNDbse;
        ATMOmax = -ATMOmin;
        ###########################################################################
        #Determine the Relative Humidity Range
        ATMOrnge = ATMOmax - ATMOmin 
        #######################################################################
        stp = np.ceil(ATMOrnge/RNDbse)*(RNDbse/TICKint)
    else:
        ATMOmin = -4; ATMOmax = 4;
        ATMOrnge = ATMOmax - ATMOmin;
        #######################################################################
        stp = np.ceil(ATMOrnge/RNDbse)*(RNDbse/TICKint)

    ATMOtick = np.arange(ATMOmin,ATMOmax+stp,stp);    
    ###########################################################################
    return ATMOmin,ATMOmax,ATMOtick