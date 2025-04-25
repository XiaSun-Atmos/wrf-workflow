#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of these functions is to compute various moisture paramaters

###############################################################################
Created on Mon Jan 28 2023
@author: dmurray (functions from badler)
"""
##################
import numpy as np

#calculate dewpoint temperature in deg C from temp and rh
def calc_dewpoint_from_RH(T,RH):
    #T: air temperature in deg C
    #RH: relative humidity in %
    E=6.112*np.exp(17.62*T/(243.12+T))
    e=RH/100*E
    Td =  243.5/(17.67/np.log(e/6.112)-1)
    return Td


#compute water vapor pressure from mixing ratio (mr in kg/kg, p in hPa)
#e in hPa
def calc_e_from_mr(MR,P):
    e=MR*P/0.622/(1+MR/0.622)
    return e

#calculates dew point (Markowski and Richardson, Eq. 225)
#e in hPa
def calc_dewpoint_from_wvp(e):
    td = 243.5/(17.67/np.log(e/6.112)-1)
    return td

#calculate dewpoint temperature in deg C from mixing ratio and pressure
def calc_dewpoint_from_MR(MR,P):
    e = calc_e_from_mr(MR,P)
    td = calc_dewpoint_from_wvp(e)
    return td

# calculate mixing ratio from temp, rh and pressure
def calc_mixingratio(temp,rh,p):
    #temp in deg C
    # rh in %
    # p in hPa
    #compute saturation water vapor pressure
    E = 6.112*np.exp(17.62*temp/(243.12+temp))
    #compute actual water vapor pressure in hPa 
    e=E*rh/100
    #return mr in g/kg
    mr=0.622*(e/(p-e))*1000
    return mr

# calculate rh from mixing ration, pressure, temp
def calc_rh_from_MR(mr,p,temp):
    #mr: in kg/kg, p: in hPa, temp: in deg C
    #rh: in %
    #first compute water vappor pressure from mixing ratio
    e=mr*p/0.622/(1+mr/0.622)
    #second compute saturation water vapor pressure from temperature
    E = 6.112*np.exp(17.62*temp/(243.12+temp))
    rh=e/E*100
    return rh

