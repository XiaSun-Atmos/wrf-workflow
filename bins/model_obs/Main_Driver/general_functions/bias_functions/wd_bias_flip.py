#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Extract Preliminary Statistics for Model 
Validation Purposes
###############################################################################
Created on Thu Dec 19 14:47:05 2019
@author: jduncan
"""
import numpy as np
import warnings

def wd_bias(MEAS,MDL):
    WDbias = MDL - MEAS;
    ###########################################################################
    """ Identify Bias Values that Need Correction """
    with np.errstate(invalid = 'ignore'):
        WDcorr_neg = np.logical_and(abs(WDbias) > 180,WDbias < 0);
        WDcorr_pos = np.logical_and(abs(WDbias) > 180,WDbias > 0);
        #######################################################################
        WDbias[WDcorr_neg] = WDbias[WDcorr_neg] + 360;
        WDbias[WDcorr_pos] = WDbias[WDcorr_pos] - 360;
        """ Old Method Defined Below 
        #WDcorr = np.where(abs(WDbias) > 180)[0]        
        ##############################################
        # Correct Values
        for c in WDcorr:
            if WDbias[c] < 0:
                WDbias[c] = WDbias[c] + 360
            else:
                WDbias[c] = WDbias[c] - 360 """
        ##############################################
        """ Calculate First-Order Statistics """
        # Root-Mean-Square Error
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            RMSE = np.sqrt(np.nanmean(WDbias.flatten()**2))
            MUb = np.nanmean(WDbias.flatten());
        #######################################################################
        # Mean Absolute Bias
        #MAD = np.nanmean(np.abs(WDbias))
        ###################################
        # return only to two decimal points
        RMSE = np.round(RMSE,decimals = 2)
        MUb = np.round(MUb,decimals = 2)
    #######################################
    return WDbias, RMSE, MUb
    
    
    
    
