#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
###############################################################################
This Function Will Output Appropriate Limits Based on the SNR Spatial Data
- The SNR data was appropriately converted to dB values from its unitless base

"""
import numpy as np
###############################################################################

def snr_sptl_ylim(SNRin,RNDbse,TICKint):
    SNRarr = [];
    ###########################################################################
    for r in SNRin:
        SNRarr = np.concatenate((SNRarr,np.ravel(r)))
    ###########################################################################
    # Do Not Establish Lower Bound of SNR Array Because of dB Units
    # SNRmin = 0;
    ###########################################################################
    # Ensure Whether Several Finite Observations Exist
    if np.count_nonzero(np.nan_to_num(SNRarr)) > 1:
        mu_ref = np.nanmean(SNRarr)
        sig_ref = np.nanstd(SNRarr)
        #######################################################################
        SNRmin = np.floor((mu_ref-2.5*sig_ref)/RNDbse)*RNDbse;
        SNRmax = np.floor((mu_ref+2.5*sig_ref)/RNDbse)*RNDbse;
        # floor used for both to highlight intricacies 
        #######################################################################
    	# if SNRmin == SNRmax -- cannot proceed so make necessary modifications
        if SNRmax == SNRmin:
    	    # modify SNRmax to allow further processing #
    	    SNRmax = np.ceil((mu_ref+2.5*sig_ref)/RNDbse)*RNDbse;
    	#######################################################################
        #SNRmin = np.floor(np.nanmin(SNRarr)/RNDbse)*RNDbse;
        #SNRmax = np.ceil(np.nanmax(SNRarr)/RNDbse)*RNDbse;
        #######################################################################
        #Determine the Relative Humidity Range
        SNRrnge = SNRmax - SNRmin 
        #######################################################################
        stp = np.ceil(SNRrnge/RNDbse)*(RNDbse/TICKint)
    else:
        SNRmin = -40; SNRmax = 10;
        SNRrnge = SNRmax - SNRmin;
        #######################################################################
        stp = np.ceil(SNRrnge/RNDbse)*(RNDbse/TICKint)

    ATMOtick = np.arange(SNRmin,SNRmax+stp,stp);    
    ###########################################################################
    return SNRmin,SNRmax,ATMOtick
