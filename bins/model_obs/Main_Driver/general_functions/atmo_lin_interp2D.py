#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Perform a Linear Interpolation of the 
Modeled Atmo Data to the Appropriate Measurement Heights
###############################################################################
Created on Tue Jan 28 14:46:04 2020
@author: dmurray
"""

import numpy as np

def atmo_lininterp2_lev2D(MDLatmo,MDLhgts,HGTSoi):
    """ Dimensions of Model Input #############################################
    ---> (hgt, time) ... (row,column) 
    ####################################################################### """
    # Perform Linear Interoplation on Atmospheric Input at each Time Stamp
    if MDLatmo.shape[1] > 0:
        ATMOint = np.asarray(np.ones((HGTSoi.shape[0],MDLatmo.shape[1]))*np.nan,dtype=float)
        #######################################################################
        for t in range(0,MDLatmo.shape[1]):
            ATMOint[:,t] = np.interp(HGTSoi,MDLhgts,MDLatmo[:,t])     
    else:
        ATMOint = np.interp(HGTSoi,MDLhgts,MDLatmo[:])    
    ###########################################################################
    return ATMOint
        
        
    
    
    
