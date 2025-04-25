#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Develop a Wind Direction Mesh Array Wherein
Data about the 0|360 Axis is Visualized Using Negative Values
###############################################################################
Created on Tue Aug 18 15:52:53 2020
@author: jamesduncan1
"""

import numpy as np
import warnings

def wd2negpos(ARRin,WDin):
    WDvals = [];
    ###########################################################################
    for r in ARRin:
        WDvals = np.concatenate((WDvals,np.ravel(r)))
    ###########################################################################
    # remove NaN Values #
    WDvals = WDvals[np.isfinite(WDvals)]
    
    # Determine Whether Correction is Required #
    NEGchck = np.count_nonzero(WDvals > 270) > 0;
    POSchck = np.count_nonzero(WDvals < 90) > 0;
    ###########################################################################
    WDout = WDin;
    if np.logical_and(NEGchck,POSchck):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ###################################################################
            WDout[WDout > 270] = WDout[WDout > 270] - 360;
    ###########################################################################
    return WDout
        



