#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Extract Preliminary Statistics for Model 
Validation Purposes
###############################################################################
Created on Thu Dec 19 10:40:27 2019
@author: jduncan
"""
import numpy as np
import warnings

def atmo_bias(bias):   
    # Root-Mean-Square Error
    with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ###################################################################
            # Mean Absolute Bias
            #MAD = np.nanmean(np.abs(bias))
            ###################################################################
            RMSE = np.sqrt(np.nanmean(bias**2))
            MUb = np.nanmean(bias)
    ###########################################################################
    # return only to two decimal points
    RMSE = np.round(RMSE,decimals = 2)
    MUb = np.round(MUb,decimals = 2)
    ###################################
    return RMSE,MUb