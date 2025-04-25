#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Limits Based on the Atmospheric Spatial
Data

"""

import numpy as np
###############################################################################

def atmo_sptl_ylim(ATMOin,RNDbse,TICKint):
    ATMOarr = [];
    ###########################################################################
    for r in ATMOin:
        ATMOarr = np.concatenate((ATMOarr,np.ravel(r)))
    ###########################################################################
    # Ensure Whether Several Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(ATMOarr)) > 1:
        #######################################################################
        MU = np.nanmean(ATMOarr);
        SIG = np.nanstd(ATMOarr);
        #######################################################################
        ATMOmax = np.ceil((MU + 2.25*SIG)/RNDbse)*RNDbse;
        ATMOmin = np.floor((MU - 2.25*SIG)/RNDbse)*RNDbse;
        #######################################################################
        # Correct Min Wind Speed if Necessary
        ATMOmin = np.nanmax([0,ATMOmin])
        #######################################################################
        #Determine the Atmospheric Input Range (Typically Wind Speed)
        ATMOrnge = ATMOmax - ATMOmin 
        #######################################################################
        stp = np.ceil(ATMOrnge/RNDbse)*(RNDbse/TICKint)
    else:
        ATMOmin = 0; ATMOmax = 30;
        ATMOrnge = ATMOmax - ATMOmin;
        #######################################################################
        stp = np.ceil(ATMOrnge/RNDbse)*(RNDbse/TICKint)

    ATMOtick = np.arange(ATMOmin,ATMOmax+stp,stp);    
    ###########################################################################
    return ATMOmin,ATMOmax,ATMOtick
