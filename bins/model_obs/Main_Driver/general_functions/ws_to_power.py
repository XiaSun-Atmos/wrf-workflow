#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to calculate wind power capacity ([0,1])
from wind speed (m/s) based on Matlab function by Laura Bianco, 06/09/2023.
###############################################################################
Created on Tue Jan 28 15:27:08 2020
@author: dmurray
"""
import numpy as np

def speed2power_offshore(speed):
    # Function to convert speed (m/s) into power capacity ([0 1])
    # Author: Laura Bianco, 05/23/2024.
    
    coeff_9 = [0.000000049702181, -0.000004308225459, 0.000159578455314, -0.003300786643873, 0.041859956205538, -0.336748546205434, 1.717104276204334, -5.342265481239600, 9.225648643661168, -6.769953313819171];
    coeff_2 = [-0.000000000000001, -0.274999999999940, 7.737499999999218];
    
    power=np.full(speed.shape,np.nan)
    power[speed <= 3] = 0
    
    speedmed1 = speed[(speed > 3) & (speed <= 16)];
    power[(speed > 3) & (speed <= 16)] = coeff_9[0]*speedmed1**9 + coeff_9[1]*speedmed1**8 +coeff_9[2]*speedmed1**7 + coeff_9[3]*speedmed1**6 + coeff_9[4]*speedmed1**5 + coeff_9[5]*speedmed1**4 + coeff_9[6]*speedmed1**3 + coeff_9[7]*speedmed1**2+coeff_9[8]*speedmed1 + coeff_9[9];
    
    power[(speed > 16) & (speed <= 24.5)] = 1;
    
    speedmed2 = speed[(speed > 24.5) & (speed < 26.5)];
    power[(speed > 24.5) & (speed < 26.5)] = coeff_2[0]*speedmed2**2 + coeff_2[1]*speedmed2 + coeff_2[2];
    
    power[(speed >= 26.5) & (speed <= 28)] = 0.45;
    power[speed > 28] = 0;
    
    power[power > 1] = 1;
    power[power < 0] = 0;

    return power
    
##def speed2power_offshore_orig(speed):
##    # Function to convert speed (m/s) into power capacity ([0 1])
##    # based on Matlab function by Laura Bianco, 06/09/2023.
##
##    coeff=[-0.000000140455305, 0.000017392357527,-0.000640292827648, 0.010778583084032,\
##           -0.094716419632650, 0.456760667910376,-1.108676014204967, 1.042883682972261]
##    power=np.full(speed.shape,np.nan)
##    power[speed <= 3] = 0
##    speedmed = speed[(speed > 3) & (speed < 14)]
##    power[(speed > 3) & (speed < 14)] = coeff[0]*speedmed**7 + coeff[1]*speedmed**6 + coeff[2]*speedmed**5 + coeff[3]*speedmed**4 + \
##        coeff[4]*speedmed**3 + coeff[5]*speedmed**2+coeff[6]*speedmed + coeff[7]
##    power[(speed >= 14) & (speed <= 25)] = 1;
##    power[speed > 25] = 0
##    power[power > 1] = 1
##    power[power < 0] = 0
##
##    return power

# Get the lists of ramp up/ramp down line.  Has a problem when up/down are consecutive because plot will connect them.
def ramp_up_down(wind_speed, threshold=4):
    speed_ramp_up = np.full(len(wind_speed),np.nan);
    speed_ramp_down = np.full(len(wind_speed),np.nan);
    for i in range(len(wind_speed)-1):
        if ((wind_speed[i+1] >= wind_speed[i]) & (wind_speed[i+1]-wind_speed[i] >= threshold)):
            speed_ramp_up[i] = wind_speed[i]
            speed_ramp_up[i+1] = wind_speed[i+1]
        elif ((wind_speed[i+1] <= wind_speed[i]) & (wind_speed[i+1]-wind_speed[i] <= -threshold)):
            speed_ramp_down[i] = wind_speed[i]
            speed_ramp_down[i+1] = wind_speed[i+1]
    return speed_ramp_up, speed_ramp_down

# Get the lists of ramp up/ramp down line segments
def ramp_up_down_segs(time, wind_speed, threshold=4):
    speed_ramp_up_list = []
    speed_ramp_down_list = []
    for i in range(len(wind_speed)-1):
        if ((wind_speed[i+1] >= wind_speed[i]) & (wind_speed[i+1]-wind_speed[i] >= threshold)):
            speed_ramp_up_list.append(((time[i],wind_speed[i]),(time[i+1],wind_speed[i+1])))
        elif ((wind_speed[i+1] <= wind_speed[i]) & (wind_speed[i+1]-wind_speed[i] <= -threshold)):
            speed_ramp_down_list.append(((time[i],wind_speed[i]),(time[i+1],wind_speed[i+1])))
    return speed_ramp_up_list, speed_ramp_down_list
