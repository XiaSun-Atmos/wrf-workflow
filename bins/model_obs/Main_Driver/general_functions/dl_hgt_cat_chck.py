#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Ensure that the Doppler Lidar Array Heights
Have the Same Measurement Height Array (This is Required for Concatenation)
###############################################################################
Created on Tue Jan  5 16:31:40 2021
@author: jamesduncan1
"""

import numpy as np

def dl_hgt_chck(pre,day,pst,ext):
    """ Define Initial Conditions
    ####################################################################### """
    HGTchck = True;
    ##########################
    PREexist = len(pre) > 0;
    DAYexist = len(day) > 0;
    PSTexist = len(pst) > 0;
    EXTexist = len(ext) > 0;
    ##########################################
    # Combine all Non-Zero Length Tuple Arrays
    INtuple = ();
    if PREexist:
        INtuple = INtuple + (pre,);
    if DAYexist:
        INtuple = INtuple + (day,);
    if PSTexist:
        INtuple = INtuple + (pst,);
    if EXTexist:
        INtuple = INtuple + (ext,);
        
    """ Proceed with Analses if Tuple Size Indicates Multiple Days Considered
    ####################################################################### """
    if len(INtuple) > 1: #indicates more that heights must be checked
        if len(INtuple) == 2:
            H1 = np.array_equal(INtuple[0][2],INtuple[1][2]);
            if not H1:
                HGTchck = False;
        #######################################################################
        elif len(INtuple) == 3:
            H1 = np.array_equal(INtuple[0][2],INtuple[1][2]);
            H2 = np.array_equal(INtuple[1][2],INtuple[2][2]);
            if not H1 or not H2:
                HGTchck = False;
        #######################################################################
        elif len(INtuple) == 4:
            H1 = np.array_equal(INtuple[0][2],INtuple[1][2]);
            H2 = np.array_equal(INtuple[1][2],INtuple[2][2]);
            H3 = np.array_equal(INtuple[2][2],INtuple[3][2]);
            if not H1 or not H2 or not H3:
                HGTchck = False;

    """ Return Height Check Outcome 
    ####################################################################### """
    return HGTchck