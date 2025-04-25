#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Produce PSL Disdrometer Plots
###############################################################################
- Edits:
    - Corrections were made so that "No Data" Plots are No Longer Produced
    - Edits were made to incorporate the operational HRRR and RAP
    - Edits were made to refrain from using NCEP | ESRL prefixes
    
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
from psldisd_main_ingest import psldisd_ingest
#from esrl_main_ingest import esrl_ingest
from pslmdl_main_ingest import pslmdl_ingest
#########################################
#from psl_wind_plot import psl_Wplot
#from psl_wind_meas_plot import psl_Wmeas
#from psl_wind_model_plot import psl_Wmodel
#########################################
from psl_pr_plot import psl_PRplot
from psl_pr_meas_plot import psl_PRmeas
from psl_pr_model_plot import psl_PRmodel
#from psl_tpmx_plot import psl_TPMXplot
#from psl_tpmx_meas_plot import psl_TPMXmeas
#from psl_tpmx_model_plot import psl_TPMXmodel

def psldisd_driver(Soi,OUTdirs,anlyper,overwrite=False):

    """ Denote Analysis for this Datastream Channel is Starting
    ####################################################################### """
    print('----> PSL Disdrometer Data Analysis Starting (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')

    """ Denote Relevant Reference Directories """
    #############################################
    meas_dir = OUTdirs[0];
    model_dir = OUTdirs[1];
    img_dir = OUTdirs[2];
    ###################################################
    """ Define Model Names """
    HRRRv4_name = 'HRRR_v4'; # 'PSL_HRRR';
    ########################################
    RAPv5_name = 'RAP_v5'; # 'PSL_RAP';

    #HRRRv4_dirname = Soi+'out';
    #RAPv5_dirname = Soi+'out';
    #-> measurement and model output directories
    #HRRRv4_out = os.path.join(model_dir,HRRRv4_dirname);
    ######################################################
    #RAPv5_out = os.path.join(model_dir,RAPv5_dirname);
    
    """ Define Relevant Image (i.e. Pre-Processed) Prefix Information """
    #####################################################################
    Phrrr_v4 = 'psldisd_' + HRRRv4_name;
    #Phrrr_v4 = 'pslmet_' + HRRRv4_name;
    #######################################
    Prap_v5 = 'psldisd_' + RAPv5_name;
    #Prap_v5 = 'pslmet_' + RAPv5_name;
    ###########################################
    # datastream prefix oi
    Pmeas_oi = 'psldisd'
    
    """ Define Analyses Period """
    ###############################
    strt_date = anlyper[0]
    end_date = anlyper[1]
    
    """ Denote the Analysis Being Performed
    ####################################################################### """
    print('----> Analyzing PSL Disdrometer Data: ' + strt_date + ' -- ' + end_date + ')')
    
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
        HRRRv4_dirname = Sid+'out';
        RAPv5_dirname = Sid+'out';
        #-> measurement and model output directories
        HRRRv4_out = os.path.join(model_dir,HRRRv4_dirname);
        #####################################################
        RAPv5_out = os.path.join(model_dir,RAPv5_dirname);
    
        """ Define Image Output Directory """
        MEASout = os.path.join(meas_dir,Sid+'out')
        IMGout = os.path.join(img_dir,Sid+'out/'+Pmeas_oi+'/')
        if not os.path.exists(IMGout):
            os.makedirs(IMGout)
        
        #IMGout_path = os.listdir(IMGout)
        IMGout_path = []
        DATEdirs = glob.glob(os.path.join(IMGout,'20*'))
        for DATEdir in DATEdirs:
           IMGout_path.extend(os.listdir(DATEdir))


        """ Identify Dates Already Processed (Model by Model) """
        #######################################################################
        HRRRv4_img_fils = [fname for fname in IMGout_path if fname.startswith(Phrrr_v4)];
        #######################################################################
        RAPv5_img_fils = [fname for fname in IMGout_path if fname.startswith(Prap_v5)];

        """ Strip Date from the Image Files
        ################################################################### """
        # HRRR v4
        HRRRv4_img_descr = [imf.split('.') for imf in HRRRv4_img_fils]; HRRRv4_img_date = [];
        for imf in range(len(HRRRv4_img_descr)):
            HRRRv4_img_date.append([dt for dt in HRRRv4_img_descr[imf] if "20" in dt][0])
            
        # RAP v5
        RAPv5_img_descr = [imf.split('.') for imf in RAPv5_img_fils]
        RAPv5_img_date = [];
        for imf in range(len(RAPv5_img_descr)):
            RAPv5_img_date.append([dt for dt in RAPv5_img_descr[imf] if "20" in dt][0])
            
        #######################################################################
        # Remove Initialization File From Each of the Image File Names
        HRRRv4_img_date = [imf.split('_')[0] for imf in HRRRv4_img_date];
        #####################################################################
        RAPv5_img_date = [imf.split('_')[0] for imf in RAPv5_img_date];
        #######################################################################
        # Determine Which Dates to Process (Strings Must Match)
        HRRRv4_date2proc = np.setdiff1d(anly_dates,HRRRv4_img_date);
        ################################################################
        RAPv5_date2proc = np.setdiff1d(anly_dates,RAPv5_img_date);
        
        if (overwrite):
            HRRRv4_date2proc = anly_dates
            RAPv5_date2proc = anly_dates

        """ Identify Available Measurement and Model Dastastreams """
        DISDdata = [fname for fname in os.listdir(MEASout) if fname.startswith(Pmeas_oi)];
        #SFCdata = [fname for fname in os.listdir(MEASout) if fname.startswith(Pmeas_oi)];
        #######################################################################
        #### HRRRv4_data = [fname for fname in os.listdir(HRRRv4_out)];
        ##############################################################
        #### RAPv5_data = [fname for fname in os.listdir(RAPv5_out)];

        """ ###################################################################
        Loop Through HRRR v4 Dates and Perform Analyses 
        ################################################################### """
        for doi in HRRRv4_date2proc:
            year=doi[:4]
            doi_date=date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8]));

            #######################################################################
            HRRRv4_data = [fname for fname in os.listdir(os.path.join(HRRRv4_out,year)) if HRRRv4_name in fname ];

            print('----> Analyzing PSL Disdrometer Data: HRRR v4 v ' + Sid + ' (' + doi +')')
            """ Ingest and Process Measurement Data """
            MEASloc,MEASxy,MEASprecip = psldisd_ingest(MEASout,DISDdata,doi)
            #MEASloc,MEASxy,MEASwind,MEAStrh,MEASpres,_ = psldisd_ingest(MEASout,SFCdata,doi)
            # PSLloc, PSLxy, PSLwind, PSLtemp_rh, PSLpres, PSLprecip,PSLwx (WX Not Available at Each Station)
            ###################################################################
            # Load in Data from doi+1
            doi_plus = (doi_date + timedelta(days=1)).strftime('%Y%m%d')
            _,MEASxy_plus,MEASprecip_plus = psldisd_ingest(MEASout,DISDdata,doi_plus)
            #_,MEASxy_plus,MEASwind_plus,MEAStrh_plus,MEASpres_plus,_ = psldisd_ingest(MEASout,SFCdata,doi_plus)
            ###################################################################
            # Load in Data from doi-1
            doi_minus = (doi_date - timedelta(days=1)).strftime('%Y%m%d')
            _,MEASxy_minus,MEASprecip_minus = psldisd_ingest(MEASout,DISDdata,doi_minus)
            #_,MEASxy_minus,MEASwind_minus,MEAStrh_minus,MEASpres_minus,_ = psldisd_ingest(MEASout,SFCdata,doi_minus)
            ###################################################################
            """ Incorporate doi + 2 for Extended Forecasts """
            doi_ext = (doi_date + timedelta(days=2)).strftime('%Y%m%d')
            _,MEASxy_ext,MEASprecip_ext = psldisd_ingest(MEASout,DISDdata,doi_ext)
            #_,MEASxy_ext,MEASwind_ext,MEAStrh_ext,MEASpres_ext,_ = psldisd_ingest(MEASout,SFCdata,doi_ext)
            
            ###################################################################
            #-> an empty tuple is returned if file does not exist
            ###################################################################            
            """ Ingest and Process HRRR v4 Data """
            #HRRRv4_ini, HRRRv4_xy, HRRRv4_loc, HRRRv4_tpmx, HRRRv4_wind, \
            #    HRRRv4_flux, HRRRv4_rad = psldisd_ingest(HRRRv4_out,HRRRv4_data,doi)
            HRRRv4_ini, HRRRv4_xy, HRRRv4_loc, HRRRv4_pr = pslmdl_ingest(HRRRv4_out,HRRRv4_data,doi)
            ################################################################
            # Functions Return a Tuple Containing Different Initalization
            # Times (Should be 24 but this Value Varies)
            ###############################################################
            FILout = (IMGout,HRRRv4_name,Sid,doi)
            #######################################

            """ Produce PSL Precip """
            # Define Measurement Time Input Tuple
            MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus,MEASxy_ext)
            ###################################################################
            # Define Atmospheric Variable Tuples #
            MEASpr_in = (MEASprecip_minus,MEASprecip,MEASprecip_plus,MEASprecip_ext)
            #MEASw_in = (MEASwind_minus,MEASwind,MEASwind_plus,MEASwind_ext)
            #MEAStrh_in = (MEAStrh_minus,MEAStrh,MEAStrh_plus,MEAStrh_ext)
            #MEASpres_in = (MEASpres_minus,MEASpres,MEASpres_plus,MEASpres_ext)
            # Declare Measurement Input Plotting Tuples #
            MEASpr_plt = (MEASloc,MEAStme,MEASpr_in)            
            #MEASw_plt = (MEASloc,MEAStme,MEASw_in)            
            #MEAStrp_plt = (MEASloc,MEAStme,MEAStrh_in,MEASpres_in)            
            ###################################################################
            # Declare Model Input Plotting Tuples #
            #HRRRpr_plt = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_tpmx)
            HRRRpr_plt = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_pr)
            #HRRRw_plt = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_wind)
            #HRRRtpmx_plt = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_tpmx)
            ###################################################################
            # Three Plotting Options:
            #   1) Model and Measurement Data are Defined 
            #   2) Measurement Data but No Model Data
            #   3) Model Data but No Measurement Data
            ###################################################################
            MEASchck = len(MEASprecip) > 0 #require measurements for doi defined
            MODELchck = len(HRRRv4_pr) > 0 #require model for doi defined
            #MEASchck = len(MEASwind) > 0 #require measurements for doi defined
            #MODELchck = len(HRRRv4_wind) > 0 #require model for doi defined
            ###################################################################
            # MEAS and MODEL checks do not need to be replicated for the other
            # atmospheric variables within the standard PSL Met Suite
            ###################################################################
            if MEASchck and MODELchck:
                # standard plotting protocol
                psl_PRplot(MEASpr_plt,HRRRpr_plt,FILout)
                #psl_TPMXplot(MEAStrp_plt,HRRRtpmx_plt,FILout)
                #psl_Wplot(MEASw_plt,HRRRw_plt,FILout)
            elif MEASchck and not MODELchck:
                # measurements but no model data
                psl_PRmeas(MEASpr_plt,FILout)
                #psl_Wmeas(MEASw_plt,FILout)
                #psl_TPMXmeas(MEAStrp_plt,FILout)
            elif MODELchck and not MEASchck:
                # model but no measurement data
                psl_PRmodel(HRRRpr_plt,FILout)
                #psl_Wmodel(HRRRw_plt,FILout)
                #psl_TPMXmodel(HRRRtpmx_plt,FILout)
        
        """ ###################################################################
        Loop Through RAP v5 Dates and Perform Analyses 
        ################################################################### """
        for doi in RAPv5_date2proc:
            year=doi[:4]
            doi_date=date(int(doi[0:4]),int(doi[4:6]),int(doi[6:8]));

            ##############################################################
            #RAPv5_data = [fname for fname in os.listdir(os.path.join(RAPv5_out,year))];
            RAPv5_data = [fname for fname in os.listdir(os.path.join(RAPv5_out,year)) if RAPv5_name in fname ];

            print('----> Analyzing PSL Disdrometer Data: RAP v5 v ' + Sid + ' (' + doi +')')
            """ Ingest and Process Measurement Data """
            MEASloc,MEASxy,MEASprecip = psldisd_ingest(MEASout,DISDdata,doi)
            #MEASloc,MEASxy,MEASwind,MEAStrh,MEASpres,_ = psldisd_ingest(MEASout,SFCdata,doi)
            # PSLloc, PSLxy, PSLwind, PSLtemp_rh, PSLpres, PSLprecip,PSLwx
            ###################################################################
            # Load in Data from doi+1
            doi_plus = (doi_date + timedelta(days=1)).strftime('%Y%m%d')
            _,MEASxy_plus,MEASprecip_plus = psldisd_ingest(MEASout,DISDdata,doi_plus)
            #_,MEASxy_plus,MEASwind_plus,MEAStrh_plus,MEASpres_plus,_ = psldisd_ingest(MEASout,SFCdata,doi_plus)
            ###################################################################
            # Load in Data from doi-1
            doi_minus = (doi_date - timedelta(days=1)).strftime('%Y%m%d')
            _,MEASxy_minus,MEASprecip_minus = psldisd_ingest(MEASout,DISDdata,doi_minus)
            #_,MEASxy_minus,MEASwind_minus,MEAStrh_minus,MEASpres_minus,_ = psldisd_ingest(MEASout,SFCdata,doi_minus)
            ###################################################################
            """ Incorporate doi + 2 for Extended Forecasts """
            doi_ext = (doi_date + timedelta(days=2)).strftime('%Y%m%d')
            _,MEASxy_ext,MEASprecip_ext = psldisd_ingest(MEASout,DISDdata,doi_ext)
            #_,MEASxy_ext,MEASwind_ext,MEAStrh_ext,MEASpres_ext,_ = psldisd_ingest(MEASout,SFCdata,doi_ext)
            
            ###################################################################
            #-> an empty tuple is returned if file does not exist
            ###################################################################
            """ Ingest and Process RAP v5 Data """
            #RAPv5_ini, RAPv5_xy, RAPv5_loc, RAPv5_tpmx, RAPv5_wind, \
            #    RAPv5_flux, RAPv5_rad = psldisd_ingest(RAPv5_out,RAPv5_data,doi)
            RAPv5_ini, RAPv5_xy, RAPv5_loc, RAPv5_pr = pslmdl_ingest(RAPv5_out,RAPv5_data,doi)
            ################################################################
            # Functions Return a Tuple Containing Different Initalization
            # Times (Should be 24 but this Value Varies)
            ###############################################################
            FILout = (IMGout,RAPv5_name,Sid,doi)
            #######################################

            """ Produce PSL Wind Speed and Direction, Temperature, RH, and Pressure Model Comparisons """
            # Define Measurement Time Input Tuple
            MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus,MEASxy_ext)
            ###################################################################
            # Define Atmospheric Variable Tuples #
            MEASpr_in = (MEASprecip_minus,MEASprecip,MEASprecip_plus,MEASprecip_ext)
            #MEASw_in = (MEASwind_minus,MEASwind,MEASwind_plus,MEASwind_ext)
            #MEAStrh_in = (MEAStrh_minus,MEAStrh,MEAStrh_plus,MEAStrh_ext)
            #MEASpres_in = (MEASpres_minus,MEASpres,MEASpres_plus,MEASpres_ext)
            # Declare Measurement Input Plotting Tuples #
            MEASpr_plt = (MEASloc,MEAStme,MEASpr_in)            
            #MEASw_plt = (MEASloc,MEAStme,MEASw_in)            
            #MEAStrp_plt = (MEASloc,MEAStme,MEAStrh_in,MEASpres_in)            
            ###################################################################
            # Declare Model Input Plotting Tuples #
            #RAPpr_plt = (RAPv5_ini,RAPv5_loc,RAPv5_xy,RAPv5_tpmx)
            RAPpr_plt = (RAPv5_ini,RAPv5_loc,RAPv5_xy,RAPv5_pr)
            #RAPw_plt = (RAPv5_ini,RAPv5_loc,RAPv5_xy,RAPv5_wind)
            #RAPtpmx_plt = (RAPv5_ini,RAPv5_loc,RAPv5_xy,RAPv5_tpmx)
            ###################################################################
            # Three Plotting Options:
            #   1) Model and Measurement Data are Defined 
            #   2) Measurement Data but No Model Data
            #   3) Model Data but No Measurement Data
            ###################################################################
            MEASchck = len(MEASprecip) > 0 #require measurements for doi defined
            MODELchck = len(RAPv5_pr) > 0 #require model for doi defined
            #MODELchck = False  #RAP data from GSL misssing
            ###################################################################
            # MEAS and MODEL checks do not need to be replicated for the other
            # atmospheric variables within the standard PSL Met Suite
            ###################################################################
            if MEASchck and MODELchck:
                # standard plotting protocol
                psl_PRplot(MEASpr_plt,RAPpr_plt,FILout)
                #psl_TPMXplot(MEAStrp_plt,RAPtpmx_plt,FILout)
                #psl_Wplot(MEASw_plt,RAPw_plt,FILout)
            elif MEASchck and not MODELchck:
                # measurements but no model data
                psl_PRmeas(MEASpr_plt,FILout)
                #psl_Wmeas(MEASw_plt,FILout)
                #psl_TPMXmeas(MEAStrp_plt,FILout)
            elif MODELchck and not MEASchck:
                # model but no measurement data
                psl_PRmodel(RAPpr_plt,FILout)
                #psl_Wmodel(RAPw_plt,FILout)
                #psl_TPMXmodel(RAPtpmx_plt,FILout)
            ###################################################################
            # no data situation is accounted for in psl_Wmodel
            ###################################################################
            
    
    """ Denote Analysis for this Datastream Channel is Complete
    ####################################################################### """
    print('----> PSL Disdrometer Data Analysis Complete (' + datetime.datetime.now().strftime('%Y-%b-%d -- %H:%M:%S MT') + ')')
    
