#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Readily Convert Between the Meteorological 
Wind Direction and the Corresponding U and V Wind Components
###############################################################################
Created on Wed Dec 18 16:38:01 2019
@author: jduncan
"""
import numpy as np

def wd2uv(WSin,WDatmo):
    """ Convert WDin to the Mathematical Wind Direction (i.e. Relative to True
    North) """
    WDmth = 270 - WDatmo;
    with np.errstate(invalid = 'ignore'):
        try:
            WDmth[WDmth < 0] = WDmth[WDmth < 0] + 360;
        except TypeError:
            if WDmth < 0:
                WDmth = WDmth + 360
    ###########################################################################
    """ Determine U and V Scaler Wind Components """
    Uwind = WSin * np.cos(np.radians(WDmth))
    Vwind = WSin * np.sin(np.radians(WDmth))
    ###########################################################################
    return Uwind,Vwind