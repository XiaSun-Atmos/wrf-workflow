#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Perform a Linear Interpolation of the 
Modeled Wind Direction Data to the Appropriate Measurement Height
###############################################################################
Created on Tue Jan 28 15:27:08 2020
@author: jduncan
"""

import numpy as np
############################
from met_to_uv import wd2uv
from uv_to_met import uv2wd
############################

def wd_lininterp2_lev(WDmdl,MDLhgts,HGToi):
    """ Dimensions of Model Input #############################################
    ---> (time, hgt) ... (column,row) 
    ####################################################################### """
    WSscale = np.ones(len(MDLhgts))
    # ~ establish scale value for performing u/v conversion
    
    # Perform Linear Interoplation on Wind Direction Input at each Time Stamp
    if WDmdl.shape[0] > 0:
        Uint = np.asarray(np.ones(WDmdl.shape[0])*np.nan,dtype=float)
        Vint = np.asarray(np.ones(WDmdl.shape[0])*np.nan,dtype=float)
        #######################################################################
        for t in range(0,WDmdl.shape[0]):
            # Convert Input Wind Direction to U and V Wind Components
            Umdl,Vmdl = wd2uv(WSscale,WDmdl[t,:]) 
            ###################################################################
            # Perform Linear Interpolation
            Uint[t] = np.interp(HGToi,MDLhgts,Umdl)   
            Vint[t] = np.interp(HGToi,MDLhgts,Vmdl)
        #######################################################################
        # Convert Array of U and V Values ot Atmospheric Wind Direction
        WDint = uv2wd(Uint,Vint)
    else:
        # Convert Input Wind Direction to U and V Wind Components
        Umdl,Vmdl = wd2uv(WSscale,WDmdl) 
        #######################################################################
        # Perform Linear Interpolation
        Uint = np.interp(HGToi,MDLhgts,Umdl)   
        Vint = np.interp(HGToi,MDLhgts,Vmdl)
        #######################################################################
        # Convert Array of U and V Values ot Atmospheric Wind Direction
        WDint = uv2wd(Uint,Vint) 
    ###########################################################################
    return WDint
        
        
    
    
    