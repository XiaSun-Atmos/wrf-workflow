#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Determine Which Input Variable is More 
Appropriate for Analyses
--> Examine the QC Variable Input 
--> Return the Array of Interest
###############################################################################
Created on Thu Jan  2 16:20:43 2020
@author: jduncan
"""
import numpy as np

def boom_arm(Bin,Bdefault):
    """ Unpack the Input Tuple Information """
    B1input = Bin[0];
    B2input = Bin[1];
    ###########################################################
    # Input has Shape of (Data, QC, Variable, Boom Arm Name) #
    ###########################################################
    """ Examine Relative Length of QC Variables to Determine Boom Arm of Interest """
    B1qc = B1input[1];
    B2qc = B2input[1];
    ###########################################################################
    # Determine the Number of Sufficient Measurements Made on Each Boom Arm 
    B1qc_len = len(np.where(B1qc == 0)[0])
    B2qc_len = len(np.where(B2qc == 0)[0])
    
    # Examine Relative Lengths (if the same length then resort to default)
    ###########################################################################
    if B1qc_len > B2qc_len:
        #######################
        Bvar_oi = B1input[0];
        Bqc_oi = B1input[1];
        Bname_oi = B1input[2];
        #######################
    elif B2qc_len > B1qc_len:
        Bvar_oi = B1input[0];
        Bqc_oi = B1input[1];
        Bname_oi = B1input[2];
    else:
        #######################
        Bvar_oi = Bin[Bdefault][0];
        Bqc_oi = Bin[Bdefault][1];
        Bname_oi = Bin[Bdefault][2];
    ###########################################################################
    return Bvar_oi,Bqc_oi,Bname_oi
