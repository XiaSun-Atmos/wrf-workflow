#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Code is to Correct the Time Series Tuple Arrays that were ceated
- To Develop an Equally Spaced time array

###############################################################################
Created on Thu Aug  6 09:51:34 2023
@author: dmurray
"""

import numpy as np

def atmo_time_correct(TIMEarr,ATMOarr):
    """ Identify Proper Meshgrid Spacing
    ####################################################################### """
    tdiff = np.diff(TIMEarr);
    # denote the median time delta #
    tstep_oi = np.nanmedian(tdiff);
    
    ###########################################################################
    # To Mitigate Erroneous Time Steps -- Insert Dumby Time Between Large Steps
    t2corr = np.where(tdiff > 1.25*tstep_oi)[0]
    ###########################################################################
    if len(t2corr) > 0:
        #MINtime = TIMEarr[0]; MAXtime = TIMEarr[-1];
        #TIMEarr_oi = np.arange(MINtime,MAXtime+0.1,tstep_oi)
        ##############################
        # Initiate Array #
        TIMEarr_oi = TIMEarr[0:t2corr[0] + 1]; # does not include large step
        TIMEarr_oi = np.hstack((TIMEarr_oi,TIMEarr_oi[-1] + tstep_oi))
        #######################################################################
        t = 1;
        while t < t2corr.shape[0]:
            TIMEarr_oi = np.hstack((TIMEarr_oi,TIMEarr[t2corr[t-1]+1:t2corr[t] + 1]));
            TIMEarr_oi = np.hstack((TIMEarr_oi,TIMEarr_oi[-1] + tstep_oi));
            t = t+1;
        #######################################################################
        # Finalize Array # 
        TIMEarr_oi = np.hstack((TIMEarr_oi,TIMEarr[t2corr[t-1]+1:]))    
        #######################################################################
        """ Reconstruct Data Arrays Arrays
        ################################################################### """
        ATMOarr_oi = np.ones(len(TIMEarr_oi))*np.nan;
        
        for i in range(0,len(TIMEarr_oi)):
            # Identify Relevant Time Indice from Input Array #
            INDoi = np.where(abs(TIMEarr_oi[i] - TIMEarr) < 1e-5)[0];
            if len(INDoi) == 1:
                # Properly Fill Measurement Grid
                ATMOarr_oi[i] = ATMOarr[INDoi[0]]
                # time indice is all that is needed
            ###############################################################
            del INDoi
    else:
        # No Need to Reconstruct Meshgrid Array -- Simply Return Input #
        #######################################################################
        TIMEarr_oi = TIMEarr;
        ATMOarr_oi = ATMOarr;

    """ Return Relevant Data Arrays
    ####################################################################### """
    return TIMEarr_oi, ATMOarr_oi
