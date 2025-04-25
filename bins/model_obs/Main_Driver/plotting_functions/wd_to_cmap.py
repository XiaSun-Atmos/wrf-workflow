#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Return the Colormap of Interest Based on the 
Span of Wind Directions Measured
###############################################################################
Created on Mon Aug 24 16:19:42 2020
@author: jamesduncan1
"""
###############################################################################
# Define Colormaps for Wind Direciton Visualization
from cmaps import find_cmap
from cmaps import get_var_cmap

#####################################################################
cmap_default=find_cmap('sunset')
cmap_wd=get_var_cmap('wd')
def wd_cmap(wd_min,wd_max):
    # Return the Colormap of Interest
    if wd_min == 0 and wd_max == 360:
        cmap_oi = cmap_wd
    else:
        #cmap_oi = cmap_default
        cmap_oi = cmap_wd
    ###########################################################################
    return cmap_oi

