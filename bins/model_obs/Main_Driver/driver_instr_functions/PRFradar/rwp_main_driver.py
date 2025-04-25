#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Produce Profiling Radar Wind Profiler Station
Plots
- The RASS System (for Temperature Purposes) is no Longer Operated by ARM 
###############################################################################
- Edits:
    - Corrections were made so that "No Data" Plots are No Longer Produced
    - Edits were made to incorporate the operational HRRR and RAP
    - Edits were made to refrain from using NCEP | ESRL prefixes
    - Edits were made to make it easier to incorporate new models

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
from rwp_main_ingest import rwp_ingest
from esrl_main_ingest import esrl_ingest
#########################################
from rwp_wind_plot import rwp_Wplot
from rwp_wind_meas_plot import rwp_Wmeas
from rwp_wind_model_plot import rwp_Wmodel
from rwp_wind_ts_plot import rwp_WTSplot
from rwp_wind_ts_meas_plot import rwp_WTSmeas
from rwp_wind_ts_model_plot import rwp_WTSmodel
#########################################
from levs2concatenate import lev_cat_chck
import station_param_list as station_prm
import forecast_param_list as fxcst_param

def rwp_driver(Soi,OUTdirs,anlyper,frequency=915,overwrite=False):

    """ Denote Analysis for this Datastream Channel is Starting
    ####################################################################### """
    print('----> '+str(frequency)+'MHz Doppler Radar Data Analysis Starting (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')

    """ Denote Relevant Reference Directories """
    #############################################
    meas_dir = OUTdirs[0];
    model_dir = OUTdirs[1];
    img_dir = OUTdirs[2];


    print(meas_dir)
    print(model_dir)
    print(img_dir)

    ###################################################
    """ Define Model Names and Locations """
    MODEL_names = getattr(fxcst_param,'model_names')
    MODEL_dirs = getattr(fxcst_param,'model_dirs')
    ###################################################
    freq_dir='rwpwinds'+str(frequency)
    freq_pfx=freq_dir+'_'
    
    # datastream prefix oi
    Pmeas_oi = 'rwp'+str(frequency)+'windhourly'
    
    """ Define Analyses Period """
    ###############################
    strt_date = anlyper[0]
    end_date = anlyper[1]
    
    """ Denote the Analysis Being Performed
    ####################################################################### """
    print('----> Analyzing '+str(frequency)+'MHz Doppler Radar Data: ' + strt_date + ' -- ' + end_date + ')')
    
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
        MEASout = os.path.join(meas_dir,Sid+'out')
        print(Sid)
        print(meas_dir)
        print(MEASout)
        IMGout = os.path.join(img_dir,Sid+'out/'+freq_dir+'/')
        if not os.path.exists(IMGout):
            os.makedirs(IMGout)
        #######################################################################

        #IMGout_path = os.listdir(IMGout)
        IMGout_path = []
        DATEdirs = glob.glob(os.path.join(IMGout,'20*'))
        for DATEdir in DATEdirs:
           IMGout_path.extend(os.listdir(DATEdir))

        # Check to see if we need to plot the time series
        stn_params=getattr(station_prm,Sid+'params',None)
        want_ts = True
        if (stn_params):
           if 'rwp' in stn_params:
              if not ('POW' in stn_params['rwp']):
                  want_ts = False

        midx=0
        for MODEL_name in MODEL_names:
            MODEL_dirname = MODEL_dirs[midx];
            MODEL_out = os.path.join(model_dir,MODEL_dirname);
            Pmodel = freq_pfx + MODEL_name;
            print(Pmodel)

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
            PSLdata = [fname for fname in os.listdir(MEASout) if fname.startswith(Pmeas_oi)]
            print(MEASout)
            print(Pmeas_oi)
            print(PSLdata) #Xia, Oct 2024
    
            """ ###################################################################
            Loop Through Model Dates and Perform Analyses 
            ################################################################### """
            for doi in MODEL_date2proc:
                year=doi[:4]
                #######################################################################
                MODEL_data = [fname for fname in os.listdir(os.path.join(MODEL_out,year))];
    
                print('----> Analyzing '+str(frequency)+'MHz Radar Wind Profiler Data: '+ MODEL_name + ' v ' + Sid + ' (' + doi +')')
                ###################################################################
                """ Ingest and Process Measurement Data """
                RWPloc,RWPxy,_,RWPwind_low,RWPwind_hgh,lev_doi = rwp_ingest(MEASout,PSLdata,doi)
                # RWPloc, RWPxy, RWPdata, RWPlow,RWPhgh
                ###################################################################
                # Load in Data from doi+1
                doi_plus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=1)).strftime('%Y%m%d')
                _,RWPxy_plus,_,RWPwind_low_plus,RWPwind_hgh_plus,lev_plus = rwp_ingest(MEASout,PSLdata,doi_plus)
                ###################################################################
                # Load in Data from doi-1
                doi_minus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) - timedelta(days=1)).strftime('%Y%m%d')
                _,RWPxy_minus,_,RWPwind_low_minus,RWPwind_hgh_minus,lev_minus = rwp_ingest(MEASout,PSLdata,doi_minus)
                ###################################################################
                """ Incorporate doi + 2 for Extended Forecasts """
                doi_ext = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=2)).strftime('%Y%m%d')
                _,RWPxy_ext,_,RWPwind_low_ext,RWPwind_hgh_ext,lev_ext = rwp_ingest(MEASout,PSLdata,doi_ext)
                ###################################################################
                # Prep for Correct Concatenation -- i.e. Ensure Equal Measurment Heights
                # lev_doi is the proper reference
                ###################################################################
                RWPxy_minus, RWPwind_low_minus, RWPwind_hgh_minus = \
                    lev_cat_chck(lev_doi,lev_minus,RWPxy_minus,RWPwind_low_minus,RWPwind_hgh_minus);
                
                RWPxy_plus, RWPwind_low_plus, RWPwind_hgh_plus = \
                    lev_cat_chck(lev_doi,lev_plus,RWPxy_plus,RWPwind_low_plus,RWPwind_hgh_plus);
                    
                RWPxy_ext, RWPwind_low_ext, RWPwind_hgh_ext = \
                    lev_cat_chck(lev_doi,lev_ext,RWPxy_ext,RWPwind_low_ext,RWPwind_hgh_ext);
                ###################################################################
                #-> an empty tuple is returned if file does not exist
                ###################################################################            
                """ Ingest and Process Model Data """
                MODEL_ini, MODEL_xy, MODEL_loc, MODEL_trp, MODEL_wind, \
                    MODEL_flux, MODEL_rad = esrl_ingest(MODEL_out,MODEL_data,doi)
                print(MODEL_data)
                print(MODEL_wind)
                ################################################################
                # Functions Return a Tuple Containing Different Initalization
                # Times (Should be 24 but this Value Varies)
                ###############################################################
                FILout = (IMGout,MODEL_name,Sid,doi)
                #######################################
    
                """ Produce Radar Wind Profiler Wind Speed and Direction Model Comparisons """
                # Define Measurement Time Input Tuple
                RWPtme = (RWPxy_minus,RWPxy,RWPxy_plus,RWPxy_ext)
                ###################################################################
                # Define Wind Speed Input Tuple
                RWPw_low = (RWPwind_low_minus,RWPwind_low,RWPwind_low_plus,RWPwind_low_ext);
                RWPw_hgh = (RWPwind_hgh_minus,RWPwind_hgh,RWPwind_hgh_plus,RWPwind_hgh_ext);
                # Declare Measurement Input Plotting Tuples 
                RWPw_plt = (RWPloc,RWPtme,RWPw_low,RWPw_hgh)            
                ###################################################################
                # Declare Model Input Plotting Tuples #
                MODELw_plt = (MODEL_ini,MODEL_loc,MODEL_xy,MODEL_wind)
                ###################################################################
                # Three Plotting Options:
                #   1) Model and Measurement Data are Defined 
                #   2) Measurement Data but No Model Data
                #   3) Model Data but No Measurement Data
                ###################################################################
                MEASchck = len(RWPwind_low) > 0 #require measurements for doi defined
                # if one mode is defined then both are defined #
                MODELchck = len(MODEL_wind) > 0 #require model for doi defined

                print(MODELchck)
                print(MEASchck)
                ###################################################################
                # MEAS and MODEL checks do not need to be replicated for the other
                # atmospheric variables within the standard ARM Met Suite
                ###################################################################
                if MEASchck and MODELchck:
                    # standard plotting protocol
                    rwp_Wplot(RWPw_plt,MODELw_plt,FILout,frequency=frequency)
                    if want_ts:
                        rwp_WTSplot(RWPw_plt,MODELw_plt,FILout,frequency=frequency)
                elif MEASchck and not MODELchck:
                    # measurements but no model data
                    rwp_Wmeas(RWPw_plt,FILout,frequency=frequency)
                    if want_ts:
                        rwp_WTSmeas(RWPw_plt,FILout,frequency=frequency)
                elif MODELchck and not MEASchck:
                    # model but no measurement data
                    rwp_Wmodel(MODELw_plt,FILout,frequency=frequency)
                    if want_ts:
                        rwp_WTSmodel(MODELw_plt,FILout,frequency=frequency)
            
            midx += 1 

    """ Denote Analysis for this Datastream Channel is Complete
    ####################################################################### """
    print('----> '+str(frequency)+'MHz Doppler Radar Data Analysis Complete (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')
