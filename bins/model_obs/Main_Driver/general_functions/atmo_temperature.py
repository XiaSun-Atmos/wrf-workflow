#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of these functions is to compute various temperature paramaters

###############################################################################
Created on Mon Jan 28 2023
@author: dmurray (functions from badler)
"""
##################
import numpy as np

def calc_theta(temp,pres):
    #function to calculate potential temperature from temperature and pressure
    #temperature in K, pressure in hPa
    theta = temp*(1000./pres)**0.2854
    return theta

def degC_to_K(temp):
    # function to convert temperature from Celsius to Kelvin
    return temp+273.15

def K_to_degC(temp):
    # function to convert temperature from Kelvin to Celsius
    return temp-273.15
