#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Code is to Correct the Spatial Tuple Arrays that were 
Created
- To Develop an Equally Spaced Meshgrid

###############################################################################
Created on Thu Aug  6 09:51:34 2020
@author: jamesduncan1
"""

import numpy as np
import math

def rwp_sptl_correct(TIMEmesh,HGTmesh,ATMOmesh):
    """ Identify Proper Meshgrid Spacing
    ####################################################################### """
    TIMEarr = TIMEmesh[0,:];
    #print(TIMEarr)
    ########################
    tdiff = np.diff(TIMEarr);
    # denote the median time delta #
    tstep_oi = np.nanmedian(tdiff);

    
    """ Extract Height and Time Information 
    ####################################################################### """
    HGTarr_oi = HGTmesh[:,0];
    ###########################################################################
    # To Mitigate Erroneous Time Steps -- Insert Dumby Time Between Large Steps
    t2corr = np.where(tdiff > 1.25*tstep_oi)[0]
    ###########################################################################
    if len(t2corr) > 0:
        MINtime = TIMEarr[0]; MAXtime = TIMEarr[-1];
        #TIMEarr_oi = np.arange(MINtime,MAXtime+0.1,tstep_oi)
        #print(TIMEarr_oi);
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
        #print(TIMEarr_oi)
        #######################################################################
        """ Reconstruct Meshgrid Arrays
        ################################################################### """
        TIMEmesh_oi, HGTmesh_oi = np.meshgrid(TIMEarr_oi,HGTarr_oi);
        ATMOmesh_oi = np.ones(TIMEmesh_oi.shape)*np.nan;
        
        for j in range(0,len(HGTarr_oi)):
            for i in range(0,len(TIMEarr_oi)):
                # Identify Relevant Time Indice from Input Array #
                INDoi = np.where(abs(TIMEarr_oi[i] - TIMEarr) < 1e-5)[0];
                if len(INDoi) == 1:
                    # Properly Fill Measurement Grid
                    ATMOmesh_oi[:,i] = ATMOmesh[:,INDoi[0]]
                    # time indice is all that is needed
                ###############################################################
                del INDoi
  
    else:
        # No Need to Reconstruct Meshgrid Array -- Simply Return Input #
        #######################################################################
        TIMEmesh_oi = TIMEmesh;
        HGTmesh_oi = HGTmesh;
        ATMOmesh_oi = ATMOmesh;

    """ Return Relevant Data Arrays
    ####################################################################### """
    return TIMEmesh_oi, HGTmesh_oi, ATMOmesh_oi
