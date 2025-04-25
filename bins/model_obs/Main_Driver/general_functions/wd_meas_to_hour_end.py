#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Purpose of this Function is to Ingest Meaurement Data and Convert it to 
Hourly Data to Enable Comparison between the Measured and Model Data
###############################################################################
Created on Wed Dec 18 14:47:46 2019
@author: jduncan

This Script Sets up the Inequality to Handle Those Measurements whose Time 
Stamp Represents the 'End' of the Averaging Interval

The Employed Inequality is:
    FCSTtime - 1800s < MEAStime <= FCSTtime + 1800s

Applicable Datastreams:
    - Surface Met Wind Direction One-Min Estimates

"""
import numpy as np
import warnings
############################
from met_to_uv import wd2uv
from uv_to_met import uv2wd
############################
HRsec = 60*60;

def wd_meas2hr_end(WDin,WDqc,Dtime,Dhr_oi,AVGmins=60,func='mean'):
    AVGsecs=AVGmins*60  # convert to seconds
    """ Define Function Input
    ###########################################################################
    WDin = Wind Direction Data In (i.e. Measured Data)
    WDqc = Optional Quality Control Information
    Dtime = Measurement Time
    Dhr_oi = Averaging Times of Interest
    ########################################################################"""
    """ Decompose Wind Direction Data into U and V Wind Components """
    WSscale = np.ones(WDin.shape);
    Uvec,Vvec = wd2uv(WSscale,WDin)
    
    #-> preallocate averaging array
    WDmu = np.ones(Dhr_oi.shape)*np.nan
    for hod in range(0,len(Dhr_oi)):
        #######################################################################
        """ Define Relevant Averaging Period"""
        #if len(WDqc)>0:
        if WDqc and len(WDqc)>0:
            AVGper = np.where(np.logical_and((np.logical_and(Dtime > Dhr_oi[hod] - AVGsecs/2,\
                 Dtime <= Dhr_oi[hod] + AVGsecs/2)),WDqc == 0))[0]
        else:
            AVGper = np.where(np.logical_and(Dtime > Dhr_oi[hod] - AVGsecs/2,Dtime <= Dhr_oi[hod] + AVGsecs/2))[0]
        #######################################################################
        """ Define the Mean U and V Wind Components """
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',category=RuntimeWarning)
            if 'median' in func:
                Umu = np.nanmedian(Uvec[AVGper])
                Vmu = np.nanmedian(Vvec[AVGper])
            else:
                Umu = np.nanmean(Uvec[AVGper])
                Vmu = np.nanmean(Vvec[AVGper])
        """ Define the Mean Wind Direction Based on Wind Components """
        WDmu[hod] = uv2wd(Umu,Vmu)
        #######################################################################
    return WDmu
