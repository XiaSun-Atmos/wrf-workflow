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

def atmo_bound_lim(ATMOin,RNDbse,MINlimit,MAXlimit,step):
    ATMOarr = [];
    ###########################################################################

    for r in ATMOin:
        ATMOarr = np.concatenate((ATMOarr,np.ravel(r)))
    ###########################################################################
    # Ensure Whether Finite Observations Exist
    print("ATMOarr below: ")
    print(ATMOarr)
    print (np.count_nonzero(np.nan_to_num(ATMOarr)))
    ATMOarr[ATMOarr>1000000]=np.nan # Xia, 2024-12-2, added this line because there are missing values for precipitation fields for column extracts for WFIP3 subhourly runs.
    print(ATMOarr)
    if np.count_nonzero(np.nan_to_num(ATMOarr)) > 0:
        print("RNDbse: "+str(RNDbse))
        print(ATMOarr)
        ATMOmin = np.floor(np.nanmin(ATMOarr)/RNDbse)*RNDbse
        ATMOmax = np.ceil(np.nanmax(ATMOarr)/RNDbse)*RNDbse
        ATMOmin = min(MINlimit,ATMOmin); ATMOmax = max(MAXlimit,ATMOmax);  #make sure at least MINlim to MAXlim
    else: #mainly for bias purposes
        ATMOmin = MINlimit;
        ATMOmax = MAXlimit;
    
    instep = step
    print("ATMOmax: "+str(ATMOmax))
    print(ATMOmin)
    print(instep)
    print((ATMOmax-ATMOmin)/step )
    while (ATMOmax-ATMOmin)/step > 8: 
       step = step + instep

    ATMOtick = np.arange(ATMOmin,ATMOmax+1,step)
    ###########################################################################
    return ATMOmin,ATMOmax,ATMOtick
