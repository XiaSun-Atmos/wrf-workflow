#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Develop and Output the Desired Doppler Lidar 
Products from the SGP Website
#############################
To Improve the 'Effective Size' of the Images, the Following Changes were Made:
    - Doppler Lidar - DL

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Aside: Differences Between Pcolor() and Pcolormesh()
    While pcolor returns a PolyCollection, pcolormesh returns a QuadMesh. 
    The latter is more specialized for the given purpose and thus is faster. 
    It should almost always be preferred. -- plt.pcolor(X,Y,DLws)
    ###########################################################################
    - plotting methods include:
    ws_img = ws_ax.pcolormesh(X,Y,DLws,vmin=Vmin,vmax=Vmax,shading='gouraud')   
    ws_img = ws_ax.imshow(DLws,extent=IMext,vmin=Vmin,vmax=Vmax,interpolation='bilinear',\
                          origin='lower',aspect='auto')                    

Using MatPlotLib's Pcolormesh to Visualize the Wind Speed Fields.
- Would like to use some type of interpolation field (i.e. Gouraud) to improve
  optics. However, when applying Gouraud (using pcolormesh) the areas 
  neighbouring NaNs appear with some weird grey structure indicating that the 
  grid box could not be properly implemented
- This describes the issue:
    https://github.com/matplotlib/matplotlib/issues/8802

Previously Attempted to Use Quiver:
###############################################################################
    #Overlay the Quiver Field
    Q = wd_ax.quiver(X[::Qres_low,::Qres_low],Y[::Qres_low,::Qres_low],DLu[::Qres_low,::Qres_low],\
                     DLv[::Qres_low,::Qres_low],units = 'xy',pivot='mid')
    #generate quiver key
    wd_ax.quiverkey(Q,0.05,1.05,10,'10 m $\mathrm{s}^{-1}$',labelpos = 'E')
###############################################################################
    More Information on the Specifications of the Quiver Plots:
        https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/quiver_demo.html#sphx-glr-gallery-images-contours-and-fields-quiver-demo-py
        https://matplotlib.org/3.1.1/gallery/images_contours_and_fields/quiver_simple_demo.html
        https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.quiver.html
        
###############################################################################
# No QC Information was Provided for these Datastreams
        
-> Updates Made to this Script Such that SNR is Only Plottted when Wind 
-> Measurements were Available        
        
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import os,time,warnings
from datetime import datetime,timedelta,timezone
#####################################
from coord_distance import dist2coord
#####################################
from atmo_meas_to_grid import atmo_meas2grid
from wd_meas_to_grid import wd_meas2grid
######################################
from time_tuple_cat import time_cat
from atmo_spatial_tuple_cat import atmo_spatial_cat
from dl_spatial_correction import dl_sptl_correct
###################################################
from atmo_ylim_plot import atmo_ylim
from atmo_spatial_ylim import atmo_sptl_ylim
from atmo_bias_spatial_ylim import atmo_sptlbias_ylim
#####################################################
from wd_spatial_ylim import wd_sptl_ylim
from wd_bias_spatial_ylim import wd_sptlbias_ylim
######################################
from wd_bias_flip import wd_bias
from cmaps import get_var_cmap
from label_funcs import gridpoint_to_label
from format_funcs import format_lat_lon_alt
####################################
#Initiate Plotting Parameter Module
import plotting_parameters as plt_param
import station_param_list as station_param
import forecast_param_list as fxcst_param

################################################################
# Define Model Interpolation Paramter
INTPoi = plt_param.INToi;
#0 indicates nearest-neighbor
#1 indicates bilinear interpolation
###############################################################################
HRsec = 60*60;
DAYsec = 24*HRsec;
MAXhgt = 200; #limit Doppler Lidar Plots to the Lowest 200 m

def dl_Wplot(MEASin,MDLin,OUTinfo):

    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth = os.path.join(OUTinfo[0], DATE);
    if not os.path.exists(IMGpth):
        os.makedirs(IMGpth)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext', DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_model = os.path.join(OUTinfo[0],'no_model', DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    
    ###########################################################################
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=image.imread(plt_param.NOAAlogo);
    FILpfx='dlwinds_'
    ATTRstr = "Credit: University of Colorado Boulder/ATOC"

    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    EXTlen = getattr(fxcst_param,MDLname.lower()+'_ext')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #    EXTlen = fxcst_param.hrrr_v4_ext;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    #    EXTlen = fxcst_param.rap_v5_ext;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    dl_Wplot(DLw_plt,HRRRw_plt,FILout)
    ###########################################################################
    MDLin ::  HRRRw_plt = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_wind)
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    ####################################################################### """
    ###########################################################################
    """ Unpack Modeled Wind Data
    ###########################
    MDLin :: 
            HRRRin = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_wind)
            RAPin = (ESRLrap_ini,ESRLrap_loc,ESRLrap_xy,ESRLrap_wind)
            #-> each embedded tuple has length of initialization files
            ###################################################################
            MDLloc = (PSLsite,PSLlon,PSLlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLwind = depends on interpolation choice (see below)
                    ESRLw_o = (WSo,WDo,WSsfc_o,WDsfc_o,Uo,Usfc_o,Vo,Vsfc_o,Wo)
                    ESRLw_int = (WSint,WDint,WDsfc_o,WDsfc_int,Uint,Usfc_int,Vint,Vsfc_int,Wint)   
                    ###########################################################################
                    ESRLw = (ESRLw_o,ESRLw_int)
    ####################################################################### """
    # Model Initialization Times #        
    MDLini = MDLin[0]; #outputs model initialization time for each ini hour
    MDLxy = MDLin[2]; #outputs forecast time and measurement height
    MDLwind = MDLin[3] #outputs wind ini hour wind data (expand in loop)
    ###########################################################################
    MDLloc = MDLin[1][0]; #just store first initialization hour
    MDLloc_names = MDLloc[0];
    MDLlon_sites = MDLloc[1]; MDLlat_sites = MDLloc[2]; # desired grid point
    ######################################################################
    MDLlon_grids = MDLloc[3]; MDLlat_grids = MDLloc[4];# model grid point referenced
    MDLelev_grids = MDLloc[5]
    MDLdist2sites = MDLloc[6]

    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
    ###########################################################################
    # measurement data not defined -- load in location from station parameter file
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat')
    MEASalt = getattr(station_param,Sid+'alt')
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    ###########################################################################
    # Get alternate grid point if available
    MEASmdlid = getattr(station_param,Sid+'mdlID',"") # mdl ID for multiple points
    ALTsite = getattr(station_param,Sid+'altSite',None)
    ALTmdlid = getattr(station_param,Sid+'altmdlID',"") # mdl ID for alternate point
    ALTmdls = getattr(station_param,Sid+'altmdls',[]) # mdl ID for alternate point
    #MEAStitle = MEASname+' ('+str(MEASalt)+' m)'
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    ###########################################################################
    # Identify the Nearest Model Grid Point to the Site Of Interest
    Grefs = []
    if (INTPoi == 0):
        Gref_meas, Gref_dist = dist2coord((MEASlat,MEASlon),(MDLlat_grids,MDLlon_grids))
    else:
        Gref_meas, Gref_dist = dist2coord((MEASlat,MEASlon),(MDLlat_sites,MDLlon_sites))
    Grefs.append(Gref_meas)
    ###########################################################################
    # If an alternate site has been requested, find the model data closest to that.
    if (ALTsite and (MDLname in ALTmdls)):
        ALTlat = ALTsite[0]; ALTlon = ALTsite[1];
        if (INTPoi == 0):
            Gref_alt, Gref_dist = dist2coord((ALTlat,ALTlon),(MDLlat_grids,MDLlon_grids))
        else:
            Gref_alt, Gref_dist = dist2coord((ALTlat,ALTlon),(MDLlat_sites,MDLlon_sites))
        Grefs.append(Gref_alt)

    ###########################################################################
    """ Unpack Measured Wind Information 
    ####################################
    MEASin :: DLw_plt = (DLloc,DLtme,DLw_in)   
    ---> DLloc = (DLlon,DLlat,DLalt)
    ---> DLtme = (DLxy_minus,DLxy,DLxy_plus)
    ---> DLw_in = (DLwind_minus,DLwind,DLwind_plus)
    ###########################################################################
    DLwind = (WS,WD,Uvel,Vvel,Wvel,SNR)
    ########################################################################"""
    # Declare Measurement Heights
    MEAShgt = MEASin[1][1][1];
    
    # Define Relevant Indices for Reference in Concatenation #
    WSind_cat = 0; WDind_cat = 1;
    Uind_cat = 2; Vind_cat = 3;
    ##########################################
    # SNRind_cat = 5; #non-range-corrected SNR
    SNRind_cat = 6; #range-corected SNR
    ###########################################################################
    # For Spatial Concatenation -- Define the Time Axis within the Input Tuple
    Tax = 1;
        
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASw_in = MEASin[2]; 
    #####################
    MEASws = atmo_spatial_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WSind_cat,Tax)
    MEASwd = atmo_spatial_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WDind_cat,Tax)
    ################################################################
    MEASu = atmo_spatial_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],Uind_cat,Tax)
    MEASv = atmo_spatial_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],Vind_cat,Tax)
    ################################################################
    
    #######################################################################
    """ Define Relevant Parameters for Doppler Lidar Plotting """
    MEAStme_base, MEAShgt_base = np.meshgrid(MEAStime,MEAShgt) #define reference meshgrid
    DLgrid = (MEAStime[0],MEAStime[-1],MEAShgt[0],MEAShgt[-1])     
    #######################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    #######################################################################
    MEAStme_mesh, MEAShgt_mesh, MEASws = dl_sptl_correct(MEAStme_base,MEAShgt_base,MEASws);
    _, _, MEASwd = dl_sptl_correct(MEAStme_base,MEAShgt_base,MEASwd);
    _, _, MEASu = dl_sptl_correct(MEAStme_base,MEAShgt_base,MEASu);
    _, _, MEASv = dl_sptl_correct(MEAStme_base,MEAShgt_base,MEASv);
            
    #######################################################################
    # set up the colormaps
    cmap_ws = get_var_cmap('ws')
    cmap_bias = get_var_cmap('bias')
    cmap_wd = get_var_cmap('wd')
    
    ###########################################################################
    # Set up some labelling for the plots based on the sites selected
    sidSfx = "";   # suffix for the station ID in the plot name
    sfxNum = 0;     # the suffix counter
    MDLstr = ""     # qualifier for the model location for plot labelling
    for Gref in Grefs:
        if (sfxNum > 0):   # the MEAS station has no suffix
           sidSfx = str(sfxNum)
        Sid = Sid+sidSfx
        #  add in the qualifier based on which point we are processing
        if (ALTsite and (MDLname in ALTmdls)): 
            if (sfxNum==0):
                MDLstr = ' ('+MEASmdlid+')'
            else:
                MDLstr = ' ('+ALTmdlid+')'

        ###########################################################################
        # Store the Reference Grid Information
        MDLlon_site =  MDLlon_sites[Gref]; MDLlat_site = MDLlat_sites[Gref]; 
        MDLlon_grid =  MDLlon_grids[Gref]; MDLlat_grid = MDLlat_grids[Gref];   
        MDLelev_grid = MDLelev_grids[Gref]
        MDLdist2site = MDLdist2sites[Gref];
        MDLloc_name = MDLloc_names[Gref]
        MDLgridinfo = gridpoint_to_label(INTPoi,sfxNum,Grefs,MDLlat_grid,MDLlon_grid,MDLelev_grid,\
                                         MDLlat_site,MDLlon_site,MEASmdlid,ALTmdlid)
        ###########################################################################
        sfxNum = sfxNum+1

        # Define the Range of Possible Initialization Times #
        INIpos = ['%02d' % it for it in range(0,24)]
        """ Initiate Measurement Model Comparison for Each Initialization Hour """
        for ini in range(0,len(INIpos)):
            # Determine Whether Model Data Exists for Said Initialization Time
            if INIpos[ini] in MDLini:
                """ Model Initialization Time Defined """
                ###################################################################
                """ Extract Relevant Model Data """
                ###########################################
                INIref = MDLini.index(INIpos[ini])
                ###########################################
                # Unpack Relevant Model Information
                MDLtime = MDLxy[INIref][0];
                ###################################################################
                # -> Convert Time (Initialization Time + Forecast Hour) to Seconds 
                # -> After Midnight
                ###################################################################
                MDLtime = int(MDLini[INIref])*HRsec + MDLtime*HRsec
                ###################################################################
    
                """ Extract Relevant Model Information """
                MDLhgt = MDLxy[INIref][1]; #model height information
                ####################################################
                MDLw_ini = MDLwind[INIref][INTPoi]
                ###################################################################
                # gridded wind information (WSo,WDo,WSsfc_o,WDsfc_o,Uo,Usfc_o,Vo,Vsfc_o,Wo)
                MDLws = MDLw_ini[0]; MDLwd = MDLw_ini[1];
                MDLu = MDLw_ini[4]; MDLv = MDLw_ini[6];
                # data is currently shaped as: (forecast_hour,sites,hgt)
                ###################################################################
                # extract from relevant grid point
                MDLws = MDLws[:,Gref,:]; MDLwd = MDLwd[:,Gref,:];
                MDLu = MDLu[:,Gref,:]; MDLv = MDLv[:,Gref,:];
                ###################################################################
                # Account for the way python interprets matrices (i.e. row x column)
                # ---> transpose inputted model data
                MDLws = np.transpose(MDLws); MDLwd = np.transpose(MDLwd);
                MDLu = np.transpose(MDLu); MDLv = np.transpose(MDLv);
                # now correctly (hgt x time) and consistent with the measured data
                ###################################################################
                """ Define Relevant Plotting Specifics 
                --> Specifics will vary depending on whether an extended forecast 
                --> is being produced
                ###################################################################
                """
                ###################################################################
                # maximum and minimum measurement height 
                HGTmin,HGTmax,HGTtick = atmo_ylim((np.asarray([0,MAXhgt])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTtick_label = [str(int(ht)) for ht in HGTtick]
                ###################################################################
                # establish barb measurement array
                BARBmeas_oi = np.zeros(MEAStme_mesh.shape,dtype = bool)
                BARBmeas_ext_oi = np.zeros(MEAStme_mesh.shape,dtype = bool)
                ###################################################################
                # define the time stamps of interest
                BARBtme_oi = HRsec*np.arange(ini + plt_param.BARBdl_x_base,(ini + FXlen),\
                              plt_param.BARBdl_x);
                BARBtme_ext_oi = HRsec*np.arange(ini + plt_param.BARBdl_x_base,(ini + EXTlen),\
                              plt_param.BARBdl_x_ext);
                ###################################################################
                # Define the Max Height Ind
                MEAShgt_max = np.where(MEAShgt > MAXhgt)[0][0] 
                MDLhgt_max = np.where(MDLhgt > MAXhgt)[0][0]
                ###################################################################
                for j in np.arange(plt_param.BARBdl_y_base,MEAShgt_max,plt_param.BARBdl_y):
                    for toi in BARBtme_oi:
                        # identify the forecast time reference indice #
                        barb_meas = np.where(np.abs(MEAStime - toi) - np.min(np.abs(MEAStime - toi)) \
                                             < 1e-5)[0][0]
                        ###########################################################
                        # only keep a single indice if equi-distant
                        ###########################################################
                        # consider if measurement delta is less than 20 min (15 min DL spacing)
                        if np.abs(MEAStime[barb_meas] - toi) < (20*60):
                            BARBmeas_oi[j,barb_meas] = True;
                    ###############################################################
                    for toi in BARBtme_ext_oi:
                        # identify the forecast time reference indice #
                        barb_meas = np.where(np.abs(MEAStime - toi) - np.min(np.abs(MEAStime - toi)) \
                                             < 1e-5)[0][0]
                        ###########################################################
                        # only keep a single indice if equi-distant
                        ###########################################################
                        # consider if measurement delta is less than 20 min (15 min DL spacing)
                        if np.abs(MEAStime[barb_meas] - toi) < (20*60):
                            BARBmeas_ext_oi[j,barb_meas] = True;
    
                ###################################################################
                if len(MDLtime) > FXlen:
                    """ Define Bounds for Extended Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLws_ext = MDLws; MDLwd_ext = MDLwd;
                    MDLu_ext = MDLu; MDLv_ext = MDLv;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_ext_mesh, MDLhgt_ext_mesh = np.meshgrid(MDLtime_ext,MDLhgt); #define reference meshgrid
                    ###############################################################
                    MEASws_ext_mu = atmo_meas2grid(MEASws,(MEAShgt_mesh,MEAStme_mesh),(MDLhgt,MDLtime_ext))
                    MEASwd_ext_mu = wd_meas2grid(MEASwd,(MEAShgt_mesh,MEAStme_mesh),(MDLhgt,MDLtime_ext))
                    ###############################################################
                    ext_meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime_ext[0], MEAStme_mesh <= MDLtime_ext[-1]))
                    ext_mdl_oi = np.logical_and(MDLhgt_ext_mesh >= MEAShgt[0], MDLhgt_ext_mesh <= MAXhgt)
                    ###############################################################
                    # pass through array indices only of interest here 
                    WSmin_ext,WSmax_ext,WStick_ext = atmo_sptl_ylim((MEASws[ext_meas_oi],\
                                    MDLws_ext[ext_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    #WDmin_ext,WDmax_ext,WDtick_ext = wd_sptl_ylim((MEASwd[ext_meas_oi],\
                    #                MDLwd_ext[ext_mdl_oi]));
                    WDmin_ext=0;WDmax_ext=360;WDtick_ext=np.arange(0,361,60);
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin_ext = MDLtime_ext[0];
                    Xmax_ext = (MDLtime_ext[0]+EXTlen*3600);
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick_ext = np.arange(MDLtime_ext[0],MDLtime_ext[0]+EXTlen*3600+1e-5,plt_param.Xdel_ext)
                    Xtick_ext_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick_ext]
                    Time_range_dates_ext = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick_ext]
                    Date_start_label_ext=Time_range_dates_ext[0].strftime('%Y-%m-%d')
                    Date_end_label_ext=Time_range_dates_ext[-1].strftime('%Y-%m-%d')
                    ###############################################################
                    ###############################################################
                    
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLws = MDLws[:,0:FXlen]; MDLwd = MDLwd[:,0:FXlen];
                    MDLu = MDLu[:,0:FXlen]; MDLv = MDLv[:,0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    ###############################################################
                    MEASws_mu = atmo_meas2grid(MEASws,(MEAShgt_mesh,MEAStme_mesh),(MDLhgt,MDLtime))
                    MEASwd_mu = wd_meas2grid(MEASwd,(MEAShgt_mesh,MEAStme_mesh),(MDLhgt,MDLtime))
                    ###############################################################
                    meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[FXlen-1]))
                    mdl_oi = np.logical_and(MDLhgt_mesh >= MEAShgt[0], MDLhgt_mesh <= MAXhgt)
                    ###############################################################
                    # pass through array indices only of interest here 
                    WSmin,WSmax,WStick = atmo_sptl_ylim((MEASws[meas_oi],MDLws[mdl_oi]),\
                                            plt_param.WSrnd_bse,plt_param.TICKint)    
                    #WDmin,WDmax,WDtick = wd_sptl_ylim((MEASwd[meas_oi],MDLwd[mdl_oi]));
                    WDmin=0;WDmax=360;WDtick=np.arange(0,361,60);
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0];
                    Xmax = MDLtime[0]+(FXlen-1)*3600;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen-1)*3600+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                    """ Define Model Barb Information 
                    ########################################################### """
                    BARBmdl_oi = np.zeros(MDLtme_mesh.shape,dtype = bool)
                    BARBmdl_ext_oi = np.zeros(MDLtme_ext_mesh.shape,dtype = bool)
                    ###############################################################
                    for j in np.arange(plt_param.BARBdl_mdl_y_base,MDLhgt_max,plt_param.BARBdl_mdl_y):
                        ###########################################################
                        # y-steps are the same whether measurement or model
                        for i in np.arange(plt_param.BARBdl_x_base,len(MDLtime),plt_param.BARBdl_x):
                            BARBmdl_oi[j,i] = True
                        ###########################################################
                        for i in np.arange(plt_param.BARBdl_x_base,len(MDLtime_ext),plt_param.BARBdl_x_ext):
                            BARBmdl_ext_oi[j,i] = True
                    
                else:
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # Know that an Extended Forecast is Non-Existent. Therefore,
                    # Simply Define meas_oi based on the Desired Forecast Time of 
                    # 86400
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    ###############################################################
                    meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[0]+(FXlen-1)*3600))
                    mdl_oi = np.logical_and(MDLhgt_mesh >= MEAShgt[0], MDLhgt_mesh <= MAXhgt)
                    ###############################################################
                    MEASws_mu = atmo_meas2grid(MEASws,(MEAShgt_mesh,MEAStme_mesh),(MDLhgt,MDLtime))
                    MEASwd_mu = wd_meas2grid(MEASwd,(MEAShgt_mesh,MEAStme_mesh),(MDLhgt,MDLtime))
                    ###############################################################
                    # pass through array indices only of interest here 
                    WSmin,WSmax,WStick = atmo_sptl_ylim((MEASws[meas_oi],MDLws[mdl_oi]),\
                                            plt_param.WSrnd_bse,plt_param.TICKint)    
                    #WDmin,WDmax,WDtick = wd_sptl_ylim((MEASwd[meas_oi],MDLwd[mdl_oi]));
                    WDmin=0;WDmax=360;WDtick=np.arange(0,361,60);
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0];
                    Xmax = MDLtime[0]+(FXlen-1)*3600;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen-1)*3600+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                    """ Define Model Barb Information 
                    ########################################################### """
                    BARBmdl_oi = np.zeros(MDLtme_mesh.shape,dtype = bool)
                    ###############################################################
                    for j in np.arange(plt_param.BARBdl_mdl_y_base,MDLhgt_max,plt_param.BARBdl_mdl_y):
                        ###########################################################
                        # y-steps are the same whether measurement or model
                        for i in np.arange(plt_param.BARBdl_x_base,len(MDLtime),plt_param.BARBdl_x):
                            BARBmdl_oi[j,i] = True
                    
                """ ###############################################################
                Wind Speed Plots (instruct to share x-axis)
                ############################################################### """
                ws_fig, (dl_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Speed Time Histories
                ##########################################
                ws_fig.suptitle('DL WS',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                #ws_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASws,vmin = WSmin,\
                ws_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASws[:-1,:-1],vmin = WSmin,\
                                            vmax = WSmax,shading='flat',cmap = cmap_ws) 
                dl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = ws_fig.colorbar(ws_dl,ax=dl_ax,ticks = WStick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                dl_ax.barbs(MEAStme_mesh[BARBmeas_oi],MEAShgt_mesh[BARBmeas_oi],MEASu[BARBmeas_oi],\
                            MEASv[BARBmeas_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                dl_ax.set_ylim(HGTmin,HGTmax); dl_ax.set_yticks(HGTtick);dl_ax.set_yticklabels(HGTtick_label);
                dl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                dl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                dl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert stn name/id as title
                T = dl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###################################################################
    
                ###################################################################
                # Plot the Model Wind Speed Time Histories
                ##########################################
                ws_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLws,vmin = WSmin,\
                                            vmax = WSmax,shading='gouraud',cmap = cmap_ws)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                mdl_ax.barbs(MDLtme_mesh[BARBmdl_oi],MDLhgt_mesh[BARBmdl_oi],MDLu[BARBmdl_oi],\
                            MDLv[BARBmdl_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WS'+MDLstr,\
                                     {'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################### """
                WSmdl_bias = MDLws - MEASws_mu;
                #######################################
                # Determine Appropriate Y-Lims for Bias
                WSb_min, WSb_max,WSb_tick = atmo_sptlbias_ylim((WSmdl_bias),plt_param.WSsptl_bias_bse,plt_param.TICKint)
                ###################################################################
                # Plot the Model Wind Speed Bias
                ##########################################
                #bias_ws = bias_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,WSmdl_bias,vmin = WSb_min,\
                bias_ws = bias_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,WSmdl_bias[:-1,:-1],vmin = WSb_min,\
                                            vmax = WSb_max,shading='flat',cmap = cmap_bias)   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = WSb_tick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ################################################################### 
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Wind Direction Plots (instruct to share x-axis)
                ############################################################### """
                wd_fig, (dl_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Direction Time Histories
                ##############################################
                wd_fig.suptitle('DL WD',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                                fontsize=plt_param.TITLEfs);
                ###################################################################
                #wd_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASwd,vmin = WDmin,\
                wd_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASwd[:-1,:-1],vmin = WDmin,\
                                            vmax = WDmax,shading='flat',cmap = cmap_wd) 
                dl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = wd_fig.colorbar(wd_dl,ax=dl_ax,ticks = WDtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                dl_ax.barbs(MEAStme_mesh[BARBmeas_oi],MEAShgt_mesh[BARBmeas_oi],MEASu[BARBmeas_oi],\
                            MEASv[BARBmeas_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                dl_ax.set_ylim(HGTmin,HGTmax); dl_ax.set_yticks(HGTtick);dl_ax.set_yticklabels(HGTtick_label);
                dl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                dl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                dl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert station info as title
                T = dl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Wind Direction Time Histories
                ##############################################
                wd_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLwd,vmin = WDmin,\
                                            vmax = WDmax,shading='gouraud',cmap = cmap_wd)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                mdl_ax.barbs(MDLtme_mesh[BARBmdl_oi],MDLhgt_mesh[BARBmdl_oi],MDLu[BARBmdl_oi],\
                            MDLv[BARBmdl_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WD'+MDLstr,\
                                     {'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Direction Biases
                ############################################################### """
                WDmdl_bias,_,_ = wd_bias(MEASwd_mu,MDLwd);
                WDb_min, WDb_max,WDb_tick = wd_sptlbias_ylim((WDmdl_bias))
                ###################################################################
                # Plot the Model Wind Direction Bias
                ##########################################
                #bias_wd = bias_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,WDmdl_bias,vmin = WDb_min,\
                bias_wd = bias_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,WDmdl_bias[:-1,:-1],vmin = WDb_min,\
                                            vmax = WDb_max,shading='flat',cmap = cmap_bias)    
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = WDb_tick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = wd_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all');
                
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    
                    """ Extended Wind Speed Plots """
                    ws_fig, (dl_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Wind Speed Time Histories
                    ##########################################
                    ws_fig.suptitle('DL WS',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #ws_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASws,vmin = WSmin_ext,\
                    ws_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASws[:-1,:-1],vmin = WSmin_ext,\
                                            vmax = WSmax_ext,shading='flat',cmap = cmap_ws) 
                    dl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = ws_fig.colorbar(ws_dl,ax=dl_ax,ticks = WStick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    # overlay barbs
                    dl_ax.barbs(MEAStme_mesh[BARBmeas_ext_oi],MEAShgt_mesh[BARBmeas_ext_oi],MEASu[BARBmeas_ext_oi],\
                                MEASv[BARBmeas_ext_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    dl_ax.set_ylim(HGTmin,HGTmax); dl_ax.set_yticks(HGTtick);dl_ax.set_yticklabels(HGTtick_label);
                    dl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    dl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                    dl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    #insert stn name/id as title
                    T = dl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Wind Speed Time Histories
                    ##########################################
                    ws_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLws_ext,\
                                               vmin = WSmin_ext,vmax = WSmax_ext,shading='gouraud',cmap = cmap_ws)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    # overlay barbs
                    mdl_ax.barbs(MDLtme_ext_mesh[BARBmdl_ext_oi],MDLhgt_ext_mesh[BARBmdl_ext_oi],MDLu_ext[BARBmdl_ext_oi],\
                                MDLv_ext[BARBmdl_ext_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' WS (m $\mathrm{s}^{-1}$)',\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Speed Biases
                    ########################################################### """
                    WSmdl_bias = MDLws_ext - MEASws_ext_mu;
                    #######################################
                    # Determine Appropriate Y-Lims for Bias
                    WSb_min, WSb_max,WSb_tick = atmo_sptlbias_ylim((WSmdl_bias),plt_param.WSsptl_bias_bse,plt_param.TICKint)
                    ###############################################################
                    # Plot the Model Wind Speed Bias
                    ##########################################
                    bias_ws = bias_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,\
                         WSmdl_bias[:-1,:-1],vmin = WSb_min,vmax = WSb_max,shading='flat',cmap = cmap_bias)   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = WSb_tick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    #bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,{'fontsize': plt_param.TICKfs});
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ############################################################### 
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###################################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                                    #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    ws_fig.savefig(FILout,dpi=300); plt.close('all');
                                         
                    """ ###########################################################
                    Wind Direction Plots (instruct to share x-axis)
                    ########################################################### """
                    wd_fig, (dl_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Wind Direction Time Histories
                    ##############################################
                    wd_fig.suptitle('DL WD',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #wd_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASwd,vmin = WDmin_ext,\
                    wd_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASwd[:-1,:-1],vmin = WDmin_ext,\
                                            vmax = WDmax_ext,shading='flat',cmap = cmap_wd); 
                    dl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = wd_fig.colorbar(wd_dl,ax=dl_ax,ticks = WDtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    # overlay barbs
                    dl_ax.barbs(MEAStme_mesh[BARBmeas_ext_oi],MEAShgt_mesh[BARBmeas_ext_oi],MEASu[BARBmeas_ext_oi],\
                                MEASv[BARBmeas_ext_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    dl_ax.set_ylim(HGTmin,HGTmax); dl_ax.set_yticks(HGTtick);dl_ax.set_yticklabels(HGTtick_label);
                    dl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    dl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                    dl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    #insert stn name/id as title
                    T = dl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Wind Direction Time Histories
                    ##############################################
                    wd_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLwd_ext,\
                                        vmin = WDmin_ext,vmax = WDmax_ext,shading='gouraud',cmap = cmap_wd);
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    # overlay barbs
                    mdl_ax.barbs(MDLtme_ext_mesh[BARBmdl_ext_oi],MDLhgt_ext_mesh[BARBmdl_ext_oi],MDLu_ext[BARBmdl_ext_oi],\
                                MDLv_ext[BARBmdl_ext_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' WD ($\mathrm{deg}^{\circ}$)',\
                                {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Direction Biases
                    ########################################################### """
                    WDmdl_bias,_,_ = wd_bias(MEASwd_ext_mu,MDLwd_ext);
                    WDb_min, WDb_max,WDb_tick = wd_sptlbias_ylim((WDmdl_bias))
                    ###############################################################
                    # Plot the Model Wind Direction Bias
                    ##########################################
                    bias_wd = bias_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,WDmdl_bias[:-1,:-1],\
                                vmin = WDb_min,vmax = WDb_max,shading='flat',cmap = cmap_bias)    
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = WDb_tick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    B = bias_ax.set_title('WD Error (Mdl - Obs)',\
                                {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###################################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                                    #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    wd_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
                    # Delete Extended Information to Preclude Incorrect Plotting
                    del MDLtime_ext
    
            else:
                """ Initialization Time Not Available """
                ###################################################################
                MDLpseudo_time = np.arange(0,(FXlen-1)*3600+1e-5,step=3600)
                MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
                MDLhgt = MEAShgt; #assumed model heights equal measurement heights
                ###################################################################
                #Define Xtick Information (Only Needs to be Defined Once)
                Xmin = MDLtime[0];
                Xmax = (MDLtime[0]+(FXlen-1)*3600);
                ###################################################################
                Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
                Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
    
                ###################################################################
                # Define Proper Indices to Determine Axes Parameters
                meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[-1]))
                ###################################################################
                # pass through array indices only of interest here 
                WSmin,WSmax,WStick = atmo_sptl_ylim((MEASws[meas_oi]),plt_param.WSrnd_bse,plt_param.TICKint)                               
                WDmin,WDmax,WDtick = wd_sptl_ylim((MEASwd[meas_oi]));
                ###################################################################
                # maximum and minimum measurement height 
                HGTmin,HGTmax,HGTtick = atmo_ylim((np.asarray([0,MAXhgt])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTtick_label = [str(int(ht)) for ht in HGTtick]
    
                # establish barb measurement array
                BARBmeas_oi = np.zeros(MEAStme_mesh.shape,dtype = bool)
                ###################################################################
                # define the time stamps of interest
                BARBtme_oi = HRsec*np.arange(ini + plt_param.BARBdl_x_base,(ini + FXlen),\
                              plt_param.BARBdl_x);
                ###################################################################
                # Define the Max Height Ind
                MEAShgt_max = np.where(MEAShgt > MAXhgt)[0][0] 
                ###################################################################
                for j in np.arange(plt_param.BARBdl_y_base,MEAShgt_max,plt_param.BARBdl_y):
                    for toi in BARBtme_oi:
                        # identify the forecast time reference indice #
                        barb_meas = np.where(np.abs(MEAStime - toi) - np.min(np.abs(MEAStime - toi)) \
                                             < 1e-5)[0][0]
                        ###########################################################
                        # only keep a single indice if equi-distant
                        ###########################################################
                        # consider if measurement delta is less than 20 min (15 min DL spacing)
                        if np.abs(MEAStime[barb_meas] - toi) < (20*60):
                            BARBmeas_oi[j,barb_meas] = True;
                
                """ ###############################################################
                Wind Speed Plots (instruct to share x-axis)
                ############################################################### """
                ws_fig, (dl_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Speed Time Histories
                ##########################################
                ws_fig.suptitle('DL WS',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                ws_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASws[:-1,:-1],vmin = WSmin,\
                                            vmax = WSmax,shading='flat', cmap = cmap_ws) 
                dl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = ws_fig.colorbar(ws_dl,ax=dl_ax,ticks = WStick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                dl_ax.barbs(MEAStme_mesh[BARBmeas_oi],MEAShgt_mesh[BARBmeas_oi],MEASu[BARBmeas_oi],\
                            MEASv[BARBmeas_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                dl_ax.set_ylim(HGTmin,HGTmax); dl_ax.set_yticks(HGTtick);dl_ax.set_yticklabels(HGTtick_label);
                dl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                dl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                dl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert stn name/id as title
                T = dl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Wind Speed Time Histories
                ##########################################
                ws_mdl = mdl_ax.imshow(MEASws*np.nan,extent=DLgrid,vmin=WSmin,vmax=WSmax,\
                                     origin='lower',aspect='auto',cmap = cmap_ws);   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WS (m $\mathrm{s}^{-1}$)',\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################### """
                # Plot the Model Wind Speed Bias
                ###################################################################
                bias_ws = bias_ax.imshow(MEASws*np.nan,extent=DLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto', cmap = cmap_bias);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Wind Direction Plots (instruct to share x-axis)
                ############################################################### """
                wd_fig, (dl_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Direction Time Histories
                ##############################################
                wd_fig.suptitle('DL WD',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                                fontsize=plt_param.TITLEfs);
                ###################################################################
                wd_dl = dl_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASwd[:-1,:-1],vmin = WDmin,\
                                            vmax = WDmax,shading='flat',cmap = cmap_wd) 
                dl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = wd_fig.colorbar(wd_dl,ax=dl_ax,ticks = WDtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                dl_ax.barbs(MEAStme_mesh[BARBmeas_oi],MEAShgt_mesh[BARBmeas_oi],MEASu[BARBmeas_oi],\
                            MEASv[BARBmeas_oi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                dl_ax.set_ylim(HGTmin,HGTmax); dl_ax.set_yticks(HGTtick);dl_ax.set_yticklabels(HGTtick_label);
                dl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                dl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                dl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert stn name/id as title
                T = dl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Wind Direction Time Histories
                ##############################################
                wd_mdl = mdl_ax.imshow(MEASwd*np.nan,extent=DLgrid,vmin=WDmin,vmax=WDmax,\
                                     origin='lower',aspect='auto',cmap = cmap_wd);    
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # Denote that Model Data was Not Available
                mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WD ($\mathrm{deg}^{\circ}$)',\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Direction Biases
                ############################################################### """
                bias_wd = bias_ax.imshow(MEASwd*np.nan,extent=DLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto', cmap = cmap_bias);    
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all');
