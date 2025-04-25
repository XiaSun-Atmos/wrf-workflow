#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""The Purpose of this Function is to Ingest Meaurement Data and Convert it to 
Hourly Data to Enable Comparison between the Measured and Model Data
###############################################################################
Created on Wed Dec 18 14:47:46 2019
@author: jduncan

This Script Sets up the Inequality to Handle Those Measurements whose Time 
Stamp Represents the 'Beginning' of the Averaging Interval

The Employed Inequality is:
    FCSTtime - 1800s <= MEAStime < FCSTtime + 1800s

Applicable Datastreams:
    - Tower Met One-Min Estimates
    - ECOR Flux Thirty-Min Estimates

"""
import numpy as np
import warnings
##################
HRsec = 60*60;

def atmo_meas2hr_beg(Din,Dqc,Dtime,Dhr_oi,AVGmins=60,func='mean'):
    AVGsecs=AVGmins*60  # convert to seconds
    """ Define Function Input
    ###########################################################################
    Din = Data In (i.e. Measured Data)
    Dqc = Optional Quality Control Information
    Dtime = Measurement Time
    Dhr_oi = Averaging Times of Interest
    ########################################################################"""
    #-> preallocate averaging array
    Dmu = np.ones(Dhr_oi.shape)*np.nan
    for hod in range(0,len(Dhr_oi)):
        #######################################################################
        """ Define Relevant Averaging Period """
        #if len(Dqc)>0:
        if Dqc and len(Dqc)>0:
            AVGper = np.where(np.logical_and((np.logical_and(Dtime >= Dhr_oi[hod] - AVGsecs/2,\
                 Dtime < Dhr_oi[hod] + AVGsecs/2)),Dqc == 0))[0]
        else:
            AVGper = np.where(np.logical_and(Dtime >= Dhr_oi[hod] - AVGsecs/2,Dtime < Dhr_oi[hod] + AVGsecs/2))[0]
        #######################################################################
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',category=RuntimeWarning)
            if 'sum' in func:
               Dmu[hod] = np.nansum(Din[AVGper])
            elif 'median' in func:
               Dmu[hod] = np.nanmedian(Din[AVGper])
            else:
               Dmu[hod] = np.nanmean(Din[AVGper])
        #######################################################################
    return Dmu
