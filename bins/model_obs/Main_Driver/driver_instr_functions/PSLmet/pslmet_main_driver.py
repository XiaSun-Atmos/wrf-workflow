#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Produce PSL SFC Meteorology Station Plots
###############################################################################
- Edits:
    - Corrections were made so that "No Data" Plots are No Longer Produced
    
###############################################################################
Created on Fri Dec 13 10:14:04 2019
@author: dmurray after jduncan
"""
import os
import glob
import numpy as np
import datetime
from datetime import timedelta, date
#########################################
from pslmet_main_ingest import pslmet_ingest
from esrl_main_ingest import esrl_ingest
#########################################
from psl_wind_plot import psl_Wplot
from psl_wind_meas_plot import psl_Wmeas
from psl_wind_model_plot import psl_Wmodel
#########################################
from psl_tpmx_plot import psl_TPMXplot
from psl_tpmx_meas_plot import psl_TPMXmeas
from psl_tpmx_model_plot import psl_TPMXmodel
#########################################
from psl_rad_plot import psl_RADplot
from psl_rad_meas_plot import psl_RADmeas
from psl_rad_model_plot import psl_RADmodel
#########################################
import station_param_list as station_prm
import forecast_param_list as fxcst_param

def pslmet_driver(Soi,OUTdirs,anlyper,overwrite=False):
    
    """ Denote Analysis for this Datastream Channel is Starting
    ####################################################################### """
    print('----> PSL MET Data Analysis Starting (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')

    """ Denote Relevant Reference Directories """
    #############################################
    meas_dir = OUTdirs[0];
    model_dir = OUTdirs[1];
    img_dir = OUTdirs[2];
    ###################################################
    """ Define Model Names and Locations """
    MODEL_names = getattr(fxcst_param,'model_names')
    MODEL_dirs = getattr(fxcst_param,'model_dirs')

    ###########################################
    # datastream prefix oi
    Pmeas_oi = 'pslmet'
    
    """ Define Analyses Period """
    ###############################
    strt_date = anlyper[0]
    end_date = anlyper[1]
    
    """ Denote the Analysis Being Performed
    ####################################################################### """
    print('----> Analyzing PSL MET Data: ' + strt_date + ' -- ' + end_date + ')')
    
    """ Develop Date Array Denoting the Relevant Analysis Period
    ####################################################################### """
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
        IMGout = os.path.join(img_dir,Sid+'out/pslmet/')
        if not os.path.exists(IMGout):
            os.makedirs(IMGout)
        
        #IMGout_path = os.listdir(IMGout)
        IMGout_path = []
        DATEdirs = glob.glob(os.path.join(IMGout,'20*'))
        for DATEdir in DATEdirs:
           IMGout_path.extend(os.listdir(DATEdir))

        midx=0
        for MODEL_name in MODEL_names:
            MODEL_dirname = MODEL_dirs[midx];
            MODEL_out = os.path.join(model_dir,MODEL_dirname);
            Pmodel = 'pslmet_' + MODEL_name;

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
            SFCdata = [fname for fname in os.listdir(MEASout) if fname.startswith(Pmeas_oi)];
            # print(MEASout)
            # print(Pmeas_oi)
            # print(SFCdata)
    
            """ ###################################################################
            Loop Through Model Dates and Perform Analyses 
            ################################################################### """
            for doi in MODEL_date2proc:
                year=doi[:4]
    
                #######################################################################
                MODEL_data = [fname for fname in os.listdir(os.path.join(MODEL_out,year))];
    
                print('----> Analyzing PSL MET Data: '+ MODEL_name + ' v ' + Sid + ' (' + doi +')')
                # print(MODEL_out)
                # print(MODEL_data)

                """ Ingest and Process Measurement Data """
                MEASloc,MEASxy,MEASwind,MEAStrh,MEASpres,MEASprecip,MEASrad = pslmet_ingest(MEASout,SFCdata,doi)
                # PSLloc, PSLxy, PSLwind, PSLtemp_rh, PSLpres, PSLprecip,PSLwx (WX Not Available at Each Station)
                ###################################################################
                # Load in Data from doi+1
                doi_plus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=1)).strftime('%Y%m%d')
                _,MEASxy_plus,MEASwind_plus,MEAStrh_plus,MEASpres_plus,MEASprecip_plus,MEASrad_plus = pslmet_ingest(MEASout,SFCdata,doi_plus)
                ###################################################################
                # Load in Data from doi-1
                doi_minus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) - timedelta(days=1)).strftime('%Y%m%d')
                _,MEASxy_minus,MEASwind_minus,MEAStrh_minus,MEASpres_minus,MEASprecip_minus,MEASrad_minus = pslmet_ingest(MEASout,SFCdata,doi_minus)
                ###################################################################
                """ Incorporate doi + 2 for Extended Forecasts """
                doi_ext = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=2)).strftime('%Y%m%d')
                _,MEASxy_ext,MEASwind_ext,MEAStrh_ext,MEASpres_ext,MEASprecip_ext,MEASrad_ext = pslmet_ingest(MEASout,SFCdata,doi_ext)

                
                ###################################################################
                #-> an empty tuple is returned if file does not exist
                ###################################################################            
                """ Ingest and Process Model Data """
                MODEL_ini, MODEL_xy, MODEL_loc, MODEL_tpmx, MODEL_wind, \
                    MODEL_flux, MODEL_rad = esrl_ingest(MODEL_out,MODEL_data,doi)
                # print(MODEL_ini)
                # print(MODEL_wind)
                ################################################################
                # Functions Return a Tuple Containing Different Initalization
                # Times (Should be 24 but this Value Varies)
                ###############################################################
                FILout = (IMGout,MODEL_name,Sid,doi)
                #######################################
    
                """ Produce PSL Wind Speed and Direction, Temperature, RH, and Pressure Model Comparisons """
                # Define Measurement Time Input Tuple
                MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus,MEASxy_ext)
                ###################################################################
                # Define Atmospheric Variable Tuples #
                MEASw_in = (MEASwind_minus,MEASwind,MEASwind_plus,MEASwind_ext)
                MEAStrh_in = (MEAStrh_minus,MEAStrh,MEAStrh_plus,MEAStrh_ext)
                MEASpres_in = (MEASpres_minus,MEASpres,MEASpres_plus,MEASpres_ext)
                MEASprecip_in = (MEASprecip_minus,MEASprecip,MEASprecip_plus,MEASprecip_ext)
                MEASrad_in = (MEASrad_minus,MEASrad,MEASrad_plus,MEASrad_ext)
                # Declare Measurement Input Plotting Tuples #
                MEASw_plt = (MEASloc,MEAStme,MEASw_in)            
                MEAStrp_plt = (MEASloc,MEAStme,MEAStrh_in,MEASpres_in,MEASprecip_in)            
                MEASrad_plt = (MEASloc,MEAStme,MEASrad_in)            
                ###################################################################
                # Declare Model Input Plotting Tuples #
                MODELw_plt = (MODEL_ini,MODEL_loc,MODEL_xy,MODEL_wind)
                MODELtpmx_plt = (MODEL_ini,MODEL_loc,MODEL_xy,MODEL_tpmx)
                MODELrad_plt = (MODEL_ini,MODEL_loc,MODEL_xy,MODEL_rad)
                ###################################################################
                # Three Plotting Options:
                #   1) Model and Measurement Data are Defined 
                #   2) Measurement Data but No Model Data
                #   3) Model Data but No Measurement Data
                ###################################################################
                MEASchck = len(MEASwind) > 0 #require measurements for doi defined
                MODELchck = len(MODEL_wind) > 0 #require model for doi defined
                print(MODELchck)
                ###################################################################
                # MEAS and MODEL checks do not need to be replicated for the other
                # atmospheric variables within the standard PSL Met Suite
                ###################################################################
                if MEASchck and MODELchck:
                    # print(FILout)
                    # standard plotting protocol
                    psl_TPMXplot(MEAStrp_plt,MODELtpmx_plt,FILout)
                    psl_Wplot(MEASw_plt,MODELw_plt,FILout)
                    psl_RADplot(MEASrad_plt,MODELrad_plt,FILout)
                elif MEASchck and not MODELchck:
                    # measurements but no model data
                    psl_TPMXmeas(MEAStrp_plt,FILout)
                    psl_Wmeas(MEASw_plt,FILout)
                    psl_RADmeas(MEASrad_plt,FILout)
                elif MODELchck and not MEASchck:
                    # model but no measurement data
                    psl_TPMXmodel(MODELtpmx_plt,FILout)
                    psl_Wmodel(MODELw_plt,FILout)
                    psl_RADmodel(MODELrad_plt,FILout)

            midx += 1

    """ Denote Analysis for this Datastream Channel is Complete
    ####################################################################### """
    print('----> PSLMET Data Analysis Complete (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')
