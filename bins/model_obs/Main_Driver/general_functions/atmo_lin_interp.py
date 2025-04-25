#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Perform a Linear Interpolation of the 
Modeled Atmo Data to the Appropriate Measurement Height
###############################################################################
Created on Tue Jan 28 14:46:04 2020
@author: jduncan
"""

import numpy as np

def atmo_lininterp2_lev(MDLatmo,MDLhgts,HGToi):
    """ Dimensions of Model Input #############################################
    ---> (time, hgt) ... (column,row) 
    ####################################################################### """
    # Perform Linear Interoplation on Atmospheric Input at each Time Stamp
    if MDLatmo.shape[0] > 0:
        ATMOint = np.asarray(np.ones(MDLatmo.shape[0])*np.nan,dtype=float)
        #######################################################################
        for t in range(0,MDLatmo.shape[0]):
            ATMOint[t] = np.interp(HGToi,MDLhgts,MDLatmo[t,:])     
    else:
        ATMOint = np.interp(HGToi,MDLhgts,MDLatmo[:])    
    ###########################################################################
    return ATMOint
        
        
    
    
    