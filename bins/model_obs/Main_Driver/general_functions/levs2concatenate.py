#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Ensure Consistency in Measurement Heights 
before Implementing the Relevant Concatenation Methods
###############################################################################
Created on Wed Aug  5 15:04:09 2020
@author: jamesduncan1
"""

def lev_cat_chck(lev2ref,lev_doi,doi_xy,doi_low_wind,doi_hgh_wind):
    # If Heights Do Not Match -- Remove Data
    if lev2ref != lev_doi:
        doi_xy = ();
        doi_low_wind = ();
        doi_hgh_wind = ();
    # Otherwise, Simply Return the Data
    ###########################################################################
    return doi_xy, doi_low_wind, doi_hgh_wind
    