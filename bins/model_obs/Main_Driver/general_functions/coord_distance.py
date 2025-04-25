#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Determine the Distance between Some Reference
Point (i.e. Lon,Lat) and an Array of Coordinates
--> The Array will Return the the Closest Two Indices and the Distance Between
--> Them

###############################################################################
Created on Wed Dec 18 09:57:46 2019
@author: jduncan
"""

from geopy import distance
import numpy as np

def dist2coord(PTref,PTanly):
    ###########################################################################
    """ Determine Distance Between Location of Interest and the Analysis Points """
    dist2pt = [distance.distance(PTref,(PTanly[0][nref],PTanly[1][nref])).m for nref in range(0,len(PTanly[0]))]

    #-> distance between points
    ###########################################################################
    """ Determine The Minimum Distance """
    ind_oi = np.where(dist2pt == np.min(dist2pt))[0][0];
    ind_dist = dist2pt[ind_oi]
    ###########################################################################
    return ind_oi,ind_dist
