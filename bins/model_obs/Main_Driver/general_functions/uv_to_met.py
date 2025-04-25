#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Readily Convert Between U and V Wind Components
and the Meteorological Wind Direction
###############################################################################
Created on Thu Dec  5 13:31:26 2019
@author: jduncan
"""

import numpy as np

def uv2wd(Uin,Vin):
    ###########################################################################
    # Convert the Input U and V Wind Vectors to the Mathematical Angle
    WDmth = np.arctan2(Vin,Uin)*180/np.pi
    ###########################################################################
    with np.errstate(invalid = 'ignore'):
        try:
            WDmth[WDmth < 0] = WDmth[WDmth < 0] +360;
        except TypeError:
            if WDmth < 0:
                WDmth = WDmth +360
    ###########################################################################
    # Convert from the Mathematical Wind Direction to the Meteorological Wind 
    # Direction
    WDatmo = 270 - WDmth;
    # Correct if Values are Less than 0
    with np.errstate(invalid = 'ignore'):
        try:
            WDatmo[WDatmo < 0] = WDatmo[WDatmo < 0] + 360;
        except TypeError:
            if WDatmo < 0:
                WDatmo = WDatmo +360
    ###########################################################################
    return WDatmo

def uv2ws(Uin,Vin):
    ###########################################################################
    # Convert the Input U and V Wind Vectors to Wind Speed
    WSatmo = np.sqrt(Uin**2 + Vin**2);
    ###########################################################################
    return WSatmo
