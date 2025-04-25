#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Produce Profiling Doppler Lidar Station Plots
###############################################################################
- Edits:
    - Edits were made to incorporate the operational HRRR and RAP
    - Edits were made to refrain from using NCEP | ESRL prefixes
    - Edits were made to make it easier to add new models
    
###############################################################################
Created on Thu Oct 24 16:37:44 2019
@author: jduncan
@modified: dmurray
"""

import os
import glob
import numpy as np
import datetime
from datetime import timedelta, date
#########################################
from dl_main_ingest import dl_ingest
from esrl_main_ingest import esrl_ingest
#########################################
from dl_hgt_cat_chck import dl_hgt_chck
#########################################
from dl_wind_plot import dl_Wplot
from dl_wind_meas_plot import dl_Wmeas
from dl_wind_model_plot import dl_Wmodel
from dl_wind_ts_plot import dl_WTSplot
from dl_wind_ts_meas_plot import dl_WTSmeas
from dl_wind_ts_model_plot import dl_WTSmodel
#########################################
import station_param_list as station_prm
import forecast_param_list as fxcst_param

def dl_driver(Soi,OUTdirs,anlyper,overwrite=False):
    """ Denote Relevant Reference Directories """
    #############################################
    meas_dir = OUTdirs[0];
    model_dir = OUTdirs[1];
    img_dir = OUTdirs[2];
    ###################################################
    """ Define Model Names and Locations """
    MODEL_names = getattr(fxcst_param,'model_names')
    MODEL_dirs = getattr(fxcst_param,'model_dirs')
    ###################################################
    dl_name='dlwinds';
    
    """ Define Analyses Period """
    ###############################
    strt_date = anlyper[0]
    end_date = anlyper[1]
    
    """ Denote the Analysis Being Performed
    ####################################################################### """
    print('----> Analyzing Doppler Lidar Data: ' + strt_date + ' -- ' + end_date + ')')
    print('----> Doppler Lidar Data Analysis Starting (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')
    
    """ Identify Relevant Dates within the Analysis Period """
    ##########################################################
    # defining date range
    anly_dates = []
    days = (date(int(end_date[0:4]),int(end_date[4:6]),int(end_date[6:8]))\
              -date(int(strt_date[0:4]),int(strt_date[4:6]),int(strt_date[6:8]))).days
    drnge = range(days+1)
    for d in drnge:
        dt = date(int(strt_date[0:4]),int(strt_date[4:6]),int(strt_date[6:8]))+timedelta(d)
        #######################################################################
        anly_dates.append(dt.strftime('%Y%m%d'))

    ###########################################################################
    """ Produce Relevant Images and Analyses for Each Station of Interest """
    for Sid in Soi:
        """ Define Image Output Directory """
        dl_Sid=getattr(station_prm,Sid+'name')  #normal station id
        dl_Pid=getattr(station_prm,Sid+'projID')

        # Check to see if we need to plot the time series
        stn_params=getattr(station_prm,Sid+'params',None)
        want_dl_ts = True
        if (stn_params):
           if 'dl' in stn_params:
              if not ('POW' in stn_params['dl']):
                  want_dl_ts = False

        MEASout = os.path.join(meas_dir,Sid+'out')
        # datastream prefix oi
        Pmeas_oi = dl_Pid.lower()+'.lidar.z03.c0.'

        IMGout = os.path.join(img_dir,Sid+'out/'+dl_name+'/')
        if not os.path.exists(IMGout):
            os.makedirs(IMGout)
        #######################################################################

        IMGout_path = []
        DATEdirs = glob.glob(os.path.join(IMGout,'20*'))
        for DATEdir in DATEdirs:
           IMGout_path.extend(os.listdir(DATEdir))

        midx=0
        for MODEL_name in MODEL_names:
            MODEL_dirname = MODEL_dirs[midx];
            MODEL_out = os.path.join(model_dir,MODEL_dirname);
            Pmodel = dl_name + '_' + MODEL_name;

            """ Identify Dates Already Processed (Model by Model) """
            #######################################################################
            MODEL_img_fils = [fname for fname in IMGout_path if fname.startswith(Pmodel)];
    
            """ Strip Date from the Image Files
            ################################################################### """
            MODEL_img_descr = [imf.split('.') for imf in MODEL_img_fils]; MODEL_img_date = [];
            for imf in range(len(MODEL_img_descr)):
                MODEL_img_date.append([dt for dt in MODEL_img_descr[imf] if "20" in dt][0])
                
            #######################################################################
            # Remove Initialization File From Each of the Image File Names
            MODEL_img_date = [imf.split('_')[0] for imf in MODEL_img_date];
            #######################################################################
            # Determine Which Dates to Process (Strings Must Match)
            MODEL_date2proc = np.setdiff1d(anly_dates,MODEL_img_date);
    
            if (overwrite):
                MODEL_date2proc = anly_dates
            
        
            """ Identify Available Measurement and Model Dastastreams """
            CUDLdata = [fname for fname in os.listdir(MEASout) if fname.startswith(Pmeas_oi)]
            #######################################################################
    
            """ ###################################################################
            Loop Through Model Dates and Perform Analyses 
            ################################################################### """
            for doi in MODEL_date2proc:
                year=doi[:4]
    
                #######################################################################
                MODEL_data = [fname for fname in os.listdir(os.path.join(MODEL_out,year))];
    
                print('----> Analyzing Doppler Lidar Data: ' + MODEL_name + ' v ' + Sid + ' (' + doi +')')
                """ Ingest and Process Measurement Data """
                DLloc,DLxy,DLwind, _, _, _ = dl_ingest(MEASout,CUDLdata,doi)
                # DLloc, DLxy, DLdata, METloc, METxy, METdata
                ###################################################################
                # Load in Data from doi+1
                doi_plus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=1)).strftime('%Y%m%d')
                _,DLxy_plus,DLwind_plus, _, _, _= dl_ingest(MEASout,CUDLdata,doi_plus)
                ###################################################################
                # Load in Data from doi-1
                doi_minus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) - timedelta(days=1)).strftime('%Y%m%d')
                _,DLxy_minus,DLwind_minus, _, _, _= dl_ingest(MEASout,CUDLdata,doi_minus)
                ###################################################################
                """ Incorporate doi + 2 for Extended Forecasts """
                doi_ext = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=2)).strftime('%Y%m%d')
                _,DLxy_ext,DLwind_ext, _, _, _= dl_ingest(MEASout,CUDLdata,doi_ext)
                
                ###################################################################
                #-> an empty tuple is returned if file does not exist
                ###################################################################            
                """ Ingest and Process Model Data """
                MODEL_ini, MODEL_xy, MODEL_loc, MODEL_trp, MODEL_wind, \
                    MODEL_flux, MODEL_rad = esrl_ingest(MODEL_out,MODEL_data,doi)
                ################################################################
                # Functions Return a Tuple Containing Different Initalization
                # Times (Should be 24 but this Value Varies)
                ###############################################################
                FILout = (IMGout,MODEL_name,Sid,doi)
                #######################################
    
                """ Produce Doppler Lidar Wind Speed and Direction Model Comparisons """
                # Define Measurement Time Input Tuple
                DLtme = (DLxy_minus,DLxy,DLxy_plus,DLxy_ext)
                ###################################################################
                # Define Wind Speed Input Tuple #
                DLw_in = (DLwind_minus,DLwind,DLwind_plus,DLwind_ext)
                # Declare Measurement Input Plotting Tuples #
                DLw_plt = (DLloc,DLtme,DLw_in)            
                ###################################################################
                # Declare Model Input Plotting Tuples #
                MODELw_plt = (MODEL_ini,MODEL_loc,MODEL_xy,MODEL_wind)
                ###################################################################
                # Three Plotting Options:
                #   1) Model and Measurement Data are Defined 
                #   2) Measurement Data but No Model Data
                #   3) Model Data but No Measurement Data
                ###################################################################
                MEASchck = len(DLwind) > 0 #require measurements for doi defined
                MODELchck = len(MODEL_wind) > 0 #require model for doi defined
                ###################################################################
                # MEAS and MODEL checks do not need to be replicated for the other
                # atmospheric variables within the standard ARM Met Suite
                ###################################################################
                """ Require Matching Measurement Height Array for All Measurement Plots """
                #HGTchck = dl_hgt_chck(DLxy_minus,DLxy,DLxy_plus,DLxy_ext);
                
                if MEASchck and MODELchck:
                    # standard plotting protocol
                    dl_Wplot(DLw_plt,MODELw_plt,FILout)
                    if (want_dl_ts):
                        dl_WTSplot(DLw_plt,MODELw_plt,FILout)
                elif MEASchck and not MODELchck:
                    # measurements but no model data
                    dl_Wmeas(DLw_plt,FILout)
                    if (want_dl_ts):
                        dl_WTSmeas(DLw_plt,FILout)
                elif MODELchck and not MEASchck:
                    # model but no measurement data
                    dl_Wmodel(MODELw_plt,FILout);
                    if (want_dl_ts):
                        dl_WTSmodel(MODELw_plt,FILout);
            midx += 1 

    """ Denote Analysis for this Datastream Channel is Complete
    ####################################################################### """
    print('----> Doppler Lidar Data Analysis Complete (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')
