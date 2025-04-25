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
    - Edits were made to make it easier to add models

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
from tropoe_main_ingest import tropoe_ingest
from esrl_main_ingest import esrl_ingest
#########################################
from tropoe_plot import tropoe_plot
from tropoe_meas_plot import tropoe_meas
from tropoe_model_plot import tropoe_model
from tropoe_ts_plot import tropoe_tsplot
from tropoe_ts_meas_plot import tropoe_tsmeas
from tropoe_ts_model_plot import tropoe_tsmodel
#########################################
from levs2concatenate import lev_cat_chck
import station_param_list as station_prm
from pslutils import get_tropoe_info
import forecast_param_list as fxcst_param

def tropoe_driver(Soi,OUTdirs,anlyper,overwrite=False,version='ASSIST'):

    """ Denote Analysis for this Datastream Channel is Starting
    ####################################################################### """
    print('----> TROPoe '+version+' Data Analysis Starting (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')

    """ Denote Relevant Reference Directories """
    #############################################
    meas_dir = OUTdirs[0];
    model_dir = OUTdirs[1];
    img_dir = OUTdirs[2];
    ###################################################
    """ Define Model Names and Locations """
    MODEL_names = getattr(fxcst_param,'model_names')
    MODEL_dirs = getattr(fxcst_param,'model_dirs')
    ########################################
    tropoe_name='tropoe'+version.lower();
    
    """ Define Analyses Period """
    ###############################
    strt_date = anlyper[0]
    end_date = anlyper[1]
    
    """ Denote the Analysis Being Performed
    ####################################################################### """
    print('----> Analyzing TROPoe '+version+' Data: (' + strt_date + ' -- ' + end_date + ')')
    
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
        # Denote the List of Data Files for this Datastream
        tropoe_dir=getattr(station_prm,Sid+version.lower());
        tropoe_Sid=getattr(station_prm,Sid+'name')  #normal station id
        tropoe_Pid=getattr(station_prm,Sid+'projID')
        MEASout,Pmeas_oi,dateindex,suffix = get_tropoe_info(station_prm.EXPname.lower(),meas_dir,\
                tropoe_dir,tropoe_Sid,tropoe_Pid,version.lower())

        """ Define Image Output Directory """
        IMGout = os.path.join(img_dir,Sid+'out/'+tropoe_name+'/')
        if not os.path.exists(IMGout):
            os.makedirs(IMGout)
        #######################################################################

        #IMGout_path = os.listdir(IMGout)
        IMGout_path = []
        DATEdirs = glob.glob(os.path.join(IMGout,'20*'))
        for DATEdir in DATEdirs:
           IMGout_path.extend(os.listdir(DATEdir))

        midx=0
        for MODEL_name in MODEL_names:
            MODEL_dirname = MODEL_dirs[midx];
            MODEL_out = os.path.join(model_dir,MODEL_dirname);
            Pmodel = tropoe_name + '_' + MODEL_name;

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
            TROPoedata = [fname for fname in os.listdir(MEASout) if ((fname.startswith(Pmeas_oi)) & (fname.endswith(suffix)))]
    
            """ ###################################################################
            Loop Through Model Dates and Perform Analyses 
            ################################################################### """
            for doi in MODEL_date2proc:
                year=doi[:4]
                #######################################################################
                MODEL_data = [fname for fname in os.listdir(os.path.join(MODEL_out,year))];
    
                print('----> Analyzing TROPoe '+version+' Data: ' + MODEL_name + ' v ' + Sid + ' (' + doi +')')
                ###################################################################
                """ Ingest and Process Measurement Data """
                TROPoeloc, TROPoexy, TROPoetemp_rh, TROPoetheta, TROPoemisc, lev_doi = tropoe_ingest(MEASout,TROPoedata,doi,suffix=suffix)
                # TROPoeloc, TROPoexy, TROPoetemp_rh, TROPoetheta
                ###################################################################
                # Load in Data from doi+1
                doi_plus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=1)).strftime('%Y%m%d')
                _, TROPoexy_plus, TROPoetemp_rh_plus, TROPoetheta_plus, TROPoemisc_plus, lev_plus = tropoe_ingest(MEASout,TROPoedata,doi_plus,suffix=suffix)
                ###################################################################
                # Load in Data from doi-1
                doi_minus = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) - timedelta(days=1)).strftime('%Y%m%d')
                _, TROPoexy_minus, TROPoetemp_rh_minus, TROPoetheta_minus, TROPoemisc_minus, lev_minus = tropoe_ingest(MEASout,TROPoedata,doi_minus,suffix=suffix)
                ###################################################################
                """ Incorporate doi + 2 for Extended Forecasts """
                doi_ext = (date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8])) + timedelta(days=2)).strftime('%Y%m%d')
                _, TROPoexy_ext, TROPoetemp_rh_ext, TROPoetheta_ext, TROPoemisc_ext, lev_ext = tropoe_ingest(MEASout,TROPoedata,doi_ext,suffix=suffix)
                ###################################################################
                # Prep for Correct Concatenation -- i.e. Ensure Equal Measurment Heights
                # lev_doi is the proper reference
                #lev_doi=0; lev_minus=0;lev_plus=0;lev_ext=0;
                ###################################################################
                #RWPxy_minus, RWPwind_low_minus, RWPwind_hgh_minus = \
                #    lev_cat_chck(lev_doi,lev_minus,RWPxy_minus,RWPwind_low_minus,RWPwind_hgh_minus);
                #
                #RWPxy_plus, RWPwind_low_plus, RWPwind_hgh_plus = \
                #    lev_cat_chck(lev_doi,lev_plus,RWPxy_plus,RWPwind_low_plus,RWPwind_hgh_plus);
                #    
                #RWPxy_ext, RWPwind_low_ext, RWPwind_hgh_ext = \
                #    lev_cat_chck(lev_doi,lev_ext,RWPxy_ext,RWPwind_low_ext,RWPwind_hgh_ext);
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
    
                """ Produce TROPoe Comparisons """
                # Define Measurement Time Input Tuple
                TROPoetme = (TROPoexy_minus,TROPoexy,TROPoexy_plus,TROPoexy_ext)
                ###################################################################
                # Define variable tuples
                TROPoetrh = (TROPoetemp_rh_minus,TROPoetemp_rh,TROPoetemp_rh_plus,TROPoetemp_rh_ext);
                TROPoetheta = (TROPoetheta_minus,TROPoetheta,TROPoetheta_plus,TROPoetheta_ext);
                TROPoemisc = (TROPoemisc_minus,TROPoemisc,TROPoemisc_plus,TROPoemisc_ext);
                # Declare Measurement Input Plotting Tuples 
                TROPoe_plt = (TROPoeloc,TROPoetme,TROPoetrh,TROPoetheta,TROPoemisc)
                ###################################################################
                # Declare Model Input Plotting Tuples #
                MODEL_plt = (MODEL_ini,MODEL_loc,MODEL_xy,MODEL_trp)
                ###################################################################
                # Three Plotting Options:
                #   1) Model and Measurement Data are Defined 
                #   2) Measurement Data but No Model Data
                #   3) Model Data but No Measurement Data
                ###################################################################
                MEASchck = len(TROPoetemp_rh) > 0 #require measurements for doi defined
                # if one mode is defined then both are defined #
                MODELchck = len(MODEL_trp) > 0 #require model for doi defined
                ###################################################################
                # MEAS and MODEL checks do not need to be replicated for the other
                # atmospheric variables within the standard ARM Met Suite
                ###################################################################
                if MEASchck and MODELchck:
                    # standard plotting protocol
                    tropoe_plot(TROPoe_plt,MODEL_plt,FILout,version=version)
                    tropoe_tsplot(TROPoe_plt,MODEL_plt,FILout,version=version)
                elif MEASchck and not MODELchck:
                    # measurements but no model data
                    tropoe_meas(TROPoe_plt,FILout,version=version)
                    tropoe_tsmeas(TROPoe_plt,FILout,version=version)
                elif MODELchck and not MEASchck:
                    # model but no measurement data
                    tropoe_model(MODEL_plt,FILout,version=version)
                    tropoe_tsmodel(MODEL_plt,FILout,version=version)
            midx += 1            
                

    """ Denote Analysis for this Datastream Channel is Complete
    ####################################################################### """
    print('----> TROPoe '+version+' Data Analysis Complete (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')
