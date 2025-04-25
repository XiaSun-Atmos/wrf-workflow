#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Where No Model Data Exists

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
#No QC Information was Provided for these Datastreams
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import os,time,warnings
from datetime import datetime,timedelta
######################################
from time_tuple_cat import time_cat
from atmo_spatial_tuple_cat import atmo_spatial_cat
from rwp_spatial_correction import rwp_sptl_correct
###################################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
from atmo_spatial_ylim import atmo_sptl_ylim
from wd_spatial_ylim import wd_sptl_ylim
from snr_spatial_ylim import snr_sptl_ylim
###################################################
from wd_to_cmap import wd_cmap
from cmaps import get_var_cmap
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

###############################################################################
# Define Desired Height Limit
MINlow = 0; MAXlow = 2500;
#MINhgh = 1250; MAXhgh = 4250;
MINhgh = 0; MAXhgh = 4250;

def rwp_Wmeas(MEASin,OUTinfo,frequency=915):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth_model = os.path.join(OUTinfo[0],'no_model',DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    ###########################################################################
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=image.imread(plt_param.NOAAlogo);
    FILpfx='rwpwinds'+str(frequency)+'_'
    
    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    rwp_Wplot(RWPw_plt,HRRRw_plt,FILout)
    ###########################################################################
    MDLin ::  HRRRw_plt = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_wind)
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    ####################################################################### """
    ###########################################################################
    """ No Modeled Wind Data to Unpack
    ###########################"""

    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure 
    ###########################################################################
    """
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat') 
    MEASalt = getattr(station_param,Sid+'alt') 
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)

    ###########################################################################
    """ Unpack Measured Wind Information 
    ####################################
    MEASin :: RWPw_plt = (RWPloc,RWPtme,RWPw_low,RWPw_hgh)  
    ---> RWPloc = (RWPlon,RWPlat,RWPalt)
    ---> RWPtme = (RWPxy_minus,RWPxy,RWPxy_plus,RWPxy_ext)
    ------> RWPw_low = (RWPwind_low_minus,RWPwind_low,RWPwind_low_plus,RWPwind_low_ext);
    ------> RWPw_hgh = (RWPwind_hgh_minus,RWPwind_hgh,RWPwind_hgh_plus,RWPwind_hgh_ext);
    ##################################################################
    # Declare Measurement Input Plotting Tuples 
    RWPw_plt = (RWPloc,RWPtme,RWPw_low,RWPw_hgh)  
    ###########################################################################
    RWPloc = (RWPlon,RWPlat,RWPalt); RWPtme = (time_meas,time_bounds,hgt);
    ###########################################################################
    RWPlow = (HGTlow,WSlow,WDlow,Uvel_low,Vvel_low,Wvel_low,SNRlow,Wlow_qc)
    RWPhgh = (HGThgh,WShgh,WDhgh,Uvel_hgh,Vvel_hgh,Wvel_hgh,SNRhgh,Whgh_qc)
    ########################################################################"""
    # Declare Measurement Heights
    MEAShgt_low = MEASin[2][1][0];
    MEAShgt_hgh = MEASin[3][1][0];
    MAXlow = MEAShgt_low[-1]
    MAXhgh = MEAShgt_hgh[-1]

    # Define Relevant Indices for Reference in Concatenation #
    WSind_cat = 1; WDind_cat = 2; 
    Uind_cat = 3; Vind_cat = 4;
    ###########################################################################
    SNRind_cat = 6; #range-corrected SNR
    Wqual_ind_cat = 7;
    ###########################################################################
    # For Spatial Concatenation -- Define the Time Axis within the Input Tuple
    Tax = 1;
        
    #######################################################################
    # set up the colormaps
    cmap_ws = get_var_cmap('ws')
    cmap_bias = get_var_cmap('bias')
    cmap_wd = get_var_cmap('wd')
    cmap_snr = get_var_cmap('snr')
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    # Define High and Low Input Tuples
    MEASw_low = MEASin[2]; 
    MEASw_hgh = MEASin[3];
    ###########################################################################
    MEASws_low = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],WSind_cat,Tax)
    MEASwd_low = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],WDind_cat,Tax)
    MEASqc_low = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],Wqual_ind_cat,Tax)
    
    MEASu_low = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],Uind_cat,Tax)
    MEASv_low = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],Vind_cat,Tax)
    
    MEASsnr_low = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],SNRind_cat,Tax)
    ###########################################################################
    MEASws_hgh = atmo_spatial_cat(MEASw_hgh[0],MEASw_hgh[1],MEASw_hgh[2],MEASw_hgh[3],WSind_cat,Tax)
    MEASwd_hgh = atmo_spatial_cat(MEASw_hgh[0],MEASw_hgh[1],MEASw_hgh[2],MEASw_hgh[3],WDind_cat,Tax)
    MEASlow_qual = atmo_spatial_cat(MEASw_low[0],MEASw_low[1],MEASw_low[2],MEASw_low[3],Wqual_ind_cat,Tax)
    
    MEASu_hgh = atmo_spatial_cat(MEASw_hgh[0],MEASw_hgh[1],MEASw_hgh[2],MEASw_hgh[3],Uind_cat,Tax)
    MEASv_hgh = atmo_spatial_cat(MEASw_hgh[0],MEASw_hgh[1],MEASw_hgh[2],MEASw_hgh[3],Vind_cat,Tax)
    
    MEASsnr_hgh = atmo_spatial_cat(MEASw_hgh[0],MEASw_hgh[1],MEASw_hgh[2],MEASw_hgh[3],SNRind_cat,Tax)
    MEAShgh_qual = atmo_spatial_cat(MEASw_hgh[0],MEASw_hgh[1],MEASw_hgh[2],MEASw_hgh[3],Wqual_ind_cat,Tax)
    ###########################################################################
    """ QC winds """
    # set wind speeds with QC > 2 to NaN
    MEASws_low[np.where(MEASlow_qual > 2)] = np.nan;
    MEASws_hgh[np.where(MEAShgh_qual > 2)] = np.nan;
    MEASwd_low[np.where(MEASlow_qual > 2)] = np.nan;
    MEASwd_hgh[np.where(MEAShgh_qual > 2)] = np.nan;
    ###########################################################################
    """ Correct SNR (i.e. only consider SNR when Wind Data was Provided/Passed
    Quality Control) """
    MEASsnr_low[np.isnan(MEASws_low)] = np.nan
    MEASsnr_hgh[np.isnan(MEASws_hgh)] = np.nan
    
    ###########################################################################
    """ Define Reference Indices for Barb Visualization """
    MEAShgt_low_max = MEAShgt_low.shape[0];
    MEAShgt_hgh_max = MEAShgt_hgh.shape[0];
    MEAShgt_low_min = np.where(MEAShgt_low >= MINlow)[0][0] 
    MEAShgt_hgh_min = np.where(MEAShgt_hgh >= MINhgh)[0][0]
    
    """ Unpack Modeled Wind Data -- No Modelled Data Provided """
    ###########################################################################    
    # Define the Range of Possible Initialization Times #
    INIpos = ['%02d' % it for it in range(0,24)]
    """ Initiate Measurement Model Comparison for Each Initialization Hour """
    for ini in range(0,len(INIpos)):
        #######################################################################
        """ Initialization Time Not Available """
        #######################################################################
        MDLpseudo_time = np.arange(0,(FXlen-1)*HRsec+1e-5,step=HRsec)
        MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
        #######################################################################
        #Define Xtick Information (Only Needs to be Defined Once)
        Xmin = MDLtime[0] - HRsec;
        Xmax = (MDLtime[0]+(FXlen-1)*HRsec) + HRsec;
        #######################################################################
        Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
        Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
        Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
        Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
        Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
        
        #######################################################################
        """ Define Relevant Parameters for Radar Wind Profiler Plotting """
        MEAStme_low_base, MEAShgt_low_base = np.meshgrid(MEAStime,MEAShgt_low) 
        MEAStme_hgh_base, MEAShgt_hgh_base = np.meshgrid(MEAStime,MEAShgt_hgh) 

        RWPgrid_low = (MEAStime[0],MEAStime[-1],MEAShgt_low[0],MEAShgt_low[-1])
        RWPgrid_hgh = (MEAStime[0],MEAStime[-1],MEAShgt_hgh[0],MEAShgt_hgh[-1])  
        # RWP Grid Used for Undefined Images 
        #######################################################################
        # Perform Correction to Create Equally Spaced Data Arrays #
        #######################################################################
        # Low Data Arrays (Only Produce Meshgrid Once)#
        MEAStme_low_mesh, MEAShgt_low_mesh, MEASws_low = rwp_sptl_correct(MEAStme_low_base,MEAShgt_low_base,MEASws_low);
        _, _, MEASwd_low = rwp_sptl_correct(MEAStme_low_base,MEAShgt_low_base,MEASwd_low);
        _, _, MEASu_low = rwp_sptl_correct(MEAStme_low_base,MEAShgt_low_base,MEASu_low);
        _, _, MEASv_low = rwp_sptl_correct(MEAStme_low_base,MEAShgt_low_base,MEASv_low);
        _, _, MEASsnr_low = rwp_sptl_correct(MEAStme_low_base,MEAShgt_low_base,MEASsnr_low);
    
        # High Data Arrays (Only Produce Meshgrid Once)#
        MEAStme_hgh_mesh, MEAShgt_hgh_mesh, MEASws_hgh = rwp_sptl_correct(MEAStme_hgh_base,MEAShgt_hgh_base,MEASws_hgh);
        _, _, MEASwd_hgh = rwp_sptl_correct(MEAStme_hgh_base,MEAShgt_hgh_base,MEASwd_hgh);
        _, _, MEASu_hgh = rwp_sptl_correct(MEAStme_hgh_base,MEAShgt_hgh_base,MEASu_hgh);
        _, _, MEASv_hgh = rwp_sptl_correct(MEAStme_hgh_base,MEAShgt_hgh_base,MEASv_hgh);
        _, _, MEASsnr_hgh = rwp_sptl_correct(MEAStme_hgh_base,MEAShgt_hgh_base,MEASsnr_hgh);
        
        #######################################################################
        # Maximum and Minimum Measurement Heights 
        HGTmin_low,HGTmax_low,HGTtick_low = atmo_ylim((np.asarray([MINlow,MAXlow])),plt_param.SHGTrnd_bse,plt_param.TICKint);
        HGTmin_hgh,HGTmax_hgh,HGTtick_hgh = atmo_ylim((np.asarray([MINhgh,MAXhgh])),plt_param.SHGTrnd_bse,plt_param.TICKint);
        HGTtick_low_label = [str(ht/1000) for ht in HGTtick_low]
        HGTtick_hgh_label = [str(ht/1000) for ht in HGTtick_hgh]
        #######################################################################
        # Define the Max Height Ind
        #MEAShgt_low_max = np.where(MEAShgt_low > MAXlow)[0][0] 
        #MEAShgt_hgh_max = np.where(MEAShgt_hgh > MAXhgh)[0][0]
        """ Establish barb array (for regular forecast lengths)
        ################################################################### """
        BARBmeas_Loi = np.zeros(MEAStme_low_mesh.shape,dtype = bool)
        BARBmeas_Hoi = np.zeros(MEAStme_hgh_mesh.shape,dtype = bool)
        #######################################################################
        # define the time stamps of interest
        BARBtme_oi = HRsec*np.arange(ini + plt_param.BARBrwp_x_base,(ini + FXlen),\
                      plt_param.BARBrwp_x);
        #######################################################################
        for toi in BARBtme_oi:
            # identify the forecast time reference indice #
            barb_meas = np.where(np.abs(MEAStime - toi) - np.min(np.abs(MEAStime - toi)) \
                                 < 1e-5)[0][0]
                ###############################################################
                # only keep a single indice if equi-distant
                ###############################################################
            # Determine Boolean for 'Low' Measurement Levels
            for j in np.arange(MEAShgt_low_min,MEAShgt_low_max,plt_param.BARBrwp_Ly):
                ###############################################################
                # consider if measurement delta is less than 6 min (1 Hr RWP spacing)
                if np.abs(MEAStime[barb_meas] - toi) < (6*60):
                    BARBmeas_Loi[j,barb_meas] = True;
            ###################################################################
            # Determine Boolean for 'High' Measurement Levels
            for j in np.arange(MEAShgt_hgh_min,MEAShgt_hgh_max,plt_param.BARBrwp_Hy):
                ###############################################################
                # consider if measurement delta is less than 6 min (1 Hr RWP spacing)
                if np.abs(MEAStime[barb_meas] - toi) < (6*60):
                    BARBmeas_Hoi[j,barb_meas] = True;      

        #######################################################################
        low_meas_oi = np.logical_and(MEAShgt_low_mesh <= MAXlow,\
                        np.logical_and(MEAStme_low_mesh >= MDLtime[0], MEAStme_low_mesh <= MDLtime[-1]));
        hgh_meas_oi = np.logical_and(MEAShgt_hgh_mesh <= MAXhgh,\
                        np.logical_and(MEAStme_hgh_mesh >= MDLtime[0], MEAStme_hgh_mesh <= MDLtime[-1]));
        #######################################################################
        # pass through array indices only of interest here 
        #WSmin_low,WSmax_low,WStick_low = atmo_sptl_ylim((MEASws_low[low_meas_oi],),plt_param.WSrnd_bse,plt_param.TICKint)    
        #WSmin_hgh,WSmax_hgh,WStick_hgh = atmo_sptl_ylim((MEASws_hgh[hgh_meas_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
        WSmin_low,WSmax_low,WStick_low = atmo_bound_lim((MEASws_low[low_meas_oi]),\
                                       plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,\
                                       plt_param.WSrnd_bse)
        WSmin_hgh,WSmax_hgh,WStick_hgh = atmo_bound_lim((MEASws_hgh[hgh_meas_oi]),\
                                       plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,\
                                       plt_param.WSrnd_bse)
        
        #WDmin_low,WDmax_low,WDtick_low = wd_sptl_ylim((MEASwd_low[low_meas_oi]));
        #WDmin_hgh,WDmax_hgh,WDtick_hgh = wd_sptl_ylim((MEASwd_hgh[hgh_meas_oi]));                                               
        WDmin_low=0;WDmax_low=360;WDtick_low=np.arange(0,361,60);
        WDmin_hgh=0;WDmax_hgh=360;WDtick_hgh=np.arange(0,361,60);
                                                                                                 
        SNRmin_low,SNRmax_low,SNRtick_low = snr_sptl_ylim((MEASsnr_low[low_meas_oi]),\
                        plt_param.SNRrnd_bse,plt_param.TICKint);
        SNRmin_hgh,SNRmax_hgh,SNRtick_hgh = snr_sptl_ylim((MEASsnr_hgh[hgh_meas_oi]),\
                        plt_param.SNRrnd_bse,plt_param.TICKint);
        
        """ Define Wind Direction Colormap Based on Array Range
        ################################################################### 
        cmap_wd_low = wd_cmap(WDmin_low,WDmax_low);
        cmap_wd_hgh = wd_cmap(WDmin_hgh,WDmax_hgh);
        """
        cmap_wd_low = cmap_wd
        cmap_wd_hgh = cmap_wd
            
        """ ###################################################################
        Wind Speed Plots (Low [i.e. High-Resolution] Mode)
        ################################################################### """
        ws_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                 height_ratios=[2,2,1],\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Lidar Wind Speed Time Histories
        ##########################################
        ws_fig.suptitle(str(frequency)+' MHz RWP WS',x = plt_param.DLmeas_x_title,\
                        y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                fontsize=plt_param.TITLEfs);
        #######################################################################
        ws_rwp = rwp_ax.pcolormesh(MEAStme_low_mesh,MEAShgt_low_mesh,MEASws_low[:-1,:-1],vmin = WSmin_low,\
                                        vmax = WSmax_low,shading='flat',cmap=cmap_ws) 
        rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
        #######################################################################
        cbar = ws_fig.colorbar(ws_rwp,ax=rwp_ax,ticks = WStick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        # overlay barbs
        rwp_ax.barbs(MEAStme_low_mesh[BARBmeas_Loi],MEAShgt_low_mesh[BARBmeas_Loi],MEASu_low[BARBmeas_Loi],\
                    MEASv_low[BARBmeas_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
        #######################################################################
        rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
        rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        # set figure y-label
        rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        #insert stn name/alt as title
        T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        ###############################

        #######################################################################
        # Plot the Model Wind Speed Time Histories
        ##########################################
        ws_mdl = mdl_ax.imshow(MEASws_low*np.nan,extent=RWPgrid_low,vmin=WSmin_low,vmax=WSmax_low,\
                             origin='lower',aspect='auto',cmap=cmap_ws); 
        mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        mdl_ax.set_ylim(HGTmin_low,HGTmax_low); mdl_ax.set_yticks(HGTtick_low);mdl_ax.set_yticklabels(HGTtick_low_label);
        mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        # Denote that Model Data was Not Available
        mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINlow,MAXlow]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        M = mdl_ax.set_title(MDLname + ' WS',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        M.set_position(plt_param.SPTLt_loc);
        
        """ ###################################################################
        Define and Plot the Model Wind Speed Biases
        ################################################################### """
        bias_ws = bias_ax.imshow(MEASws_low*np.nan,extent=RWPgrid_low,vmin = -3,vmax = 3,\
                             origin='lower',aspect='auto',cmap=cmap_bias); 
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
        #######################################################################
        bias_ax.set_ylim(HGTmin_low,HGTmax_low); bias_ax.set_yticks(HGTtick_low);bias_ax.set_yticklabels(HGTtick_low_label);
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bias_ax.set_xlim(Xmin,Xmax)
        ####################################################################### 
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        # Denote that Model Data was Not Available
        bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINlow,MAXlow]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
        ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WSlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        ws_fig.savefig(FILout,dpi=300); plt.close('all');
                             
        """ ###################################################################
        Wind Direction Plots (Low [i.e. High-Resolution] Mode)
        ################################################################### """
        wd_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                 height_ratios=[2,2,1],\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Lidar Wind Direction Time Histories
        ##############################################
        wd_fig.suptitle(str(frequency)+' MHz RWP WD',x = plt_param.DLmeas_x_title,\
                        y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
        #######################################################################
        wd_rwp = rwp_ax.pcolormesh(MEAStme_low_mesh,MEAShgt_low_mesh,MEASwd_low[:-1,:-1],vmin = WDmin_low,\
                                        vmax = WDmax_low,shading='flat',cmap=cmap_wd_low); 
        rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
        #######################################################################
        cbar = wd_fig.colorbar(wd_rwp,ax=rwp_ax,ticks = WDtick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        # overlay barbs
        rwp_ax.barbs(MEAStme_low_mesh[BARBmeas_Loi],MEAShgt_low_mesh[BARBmeas_Loi],MEASu_low[BARBmeas_Loi],\
                    MEASv_low[BARBmeas_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
        #######################################################################
        rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
        rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        # set figure y-label
        rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        #insert stn name/alt as title
        T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        ###############################

        #######################################################################
        # Plot the Model Wind Direction Time Histories
        ##############################################
        wd_mdl = mdl_ax.imshow(MEASwd_low*np.nan,extent=RWPgrid_low,vmin=WDmin_low,vmax=WDmax_low,\
                             origin='lower',aspect='auto',cmap=cmap_wd_low);  
        mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        mdl_ax.set_ylim(HGTmin_low,HGTmax_low); mdl_ax.set_yticks(HGTtick_low);mdl_ax.set_yticklabels(HGTtick_low_label);
        mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINlow,MAXlow]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        M = mdl_ax.set_title(MDLname + ' WD',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        M.set_position(plt_param.SPTLt_loc);
        
        """ ###################################################################
        Define and Plot the Model Wind Direction Biases
        ################################################################### """
        bias_wd = bias_ax.imshow(MEASwd_low*np.nan,extent=RWPgrid_low,vmin = -3,vmax = 3,\
                             origin='lower',aspect='auto',cmap=cmap_bias); 
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        bias_ax.set_ylim(HGTmin_low,HGTmax_low); bias_ax.set_yticks(HGTtick_low);bias_ax.set_yticklabels(HGTtick_low_label);
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINlow,MAXlow]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
        wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WDlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        wd_fig.savefig(FILout,dpi=300); plt.close('all');
                             
        """ ###################################################################
        Signal-to-Noise Ratio Plots (Low [i.e. High-Resolution] Mode)
        ################################################################### """
        #snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
        snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                 height_ratios=[3,1,1],\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Lidar Signal-to-Noise Time Histories
        ###############################################
        snr_fig.suptitle(str(frequency)+' MHz RWP Range-Corrected SNR',x = plt_param.DLmeas_x_title,\
                        y = plt_param.DLmeas_y_title,horizontalalignment='right',\
                fontsize=plt_param.TITLEfs);
        #######################################################################
        snr_rwp = rwp_ax.pcolormesh(MEAStme_low_mesh,MEAShgt_low_mesh,MEASsnr_low[:-1,:-1],vmin = SNRmin_low,\
                                        vmax = SNRmax_low,shading='flat',cmap=cmap_snr) 
        rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
        #######################################################################
        # overlay barbs
        rwp_ax.barbs(MEAStme_low_mesh[BARBmeas_Loi],MEAShgt_low_mesh[BARBmeas_Loi],MEASu_low[BARBmeas_Loi],\
                    MEASv_low[BARBmeas_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
        #######################################################################
        cbar = snr_fig.colorbar(snr_rwp,ax=rwp_ax,ticks = SNRtick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect3)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('SNR, dB', fontsize=plt_param.LABELfs);
        #######################################################################
        rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
        rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        # set figure y-label
        rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        #insert stn name/alt as title
        T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        ###############################          
                             
        """ ###################################################################
        Plot the Model Wind Speed Biases
        ################################################################### """
        #bias_ws = ws_ax.imshow(MEASws_low*np.nan,extent=RWPgrid_low,vmin=WSmin_low,vmax=WSmax_low,\
        bias_ws = ws_ax.imshow(MEASws_low*np.nan,extent=RWPgrid_low,vmin=-3,vmax=3,\
                             origin='lower',aspect='auto',cmap=cmap_bias); 
        ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(bias_ws,ax=ws_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
        #######################################################################
        ws_ax.set_ylim(HGTmin_low,HGTmax_low); ws_ax.set_yticks(HGTtick_low);ws_ax.set_yticklabels(HGTtick_low_label);
        ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        ws_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        ws_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINlow,MAXlow]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        B = ws_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################                     
                             
        """ ###################################################################
        Define and Plot the Model Wind Direction Biases
        ################################################################### """
        #bias_wd = wd_ax.imshow(MEASwd_low*np.nan,extent=RWPgrid_low,vmin=WSmin_low,vmax=WSmax_low,\
        bias_wd = wd_ax.imshow(MEASwd_low*np.nan,extent=RWPgrid_low,vmin=-30,vmax=30,\
                             origin='lower',aspect='auto',cmap=cmap_bias);  
        wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(bias_wd,ax=wd_ax,ticks = np.arange(-30,30.1,10),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        wd_ax.set_ylim(HGTmin_low,HGTmax_low); wd_ax.set_yticks(HGTtick_low);wd_ax.set_yticklabels(HGTtick_low_label);
        wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        wd_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        wd_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINlow,MAXlow]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        wd_ax.set_xticks(Xtick); wd_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            wd_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        wd_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        wd_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        B = wd_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
        snr_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
        #######################################################################
        wd_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        wd_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_SNRlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        snr_fig.savefig(FILout,dpi=300); plt.close('all');  
        #######################################################################
        #######################################################################
        #######################################################################
        """ ###############################################################
        Wind Speed Plots (High [i.e. Low-Resolution] Mode)
        ############################################################### """
        ws_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                 height_ratios=[2,2,1],\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Lidar Wind Speed Time Histories
        ##########################################
        ws_fig.suptitle(str(frequency)+' MHz RWP WS',x = plt_param.DLmeas_x_title,\
                        y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                fontsize=plt_param.TITLEfs);
        #######################################################################
        #ws_rwp = rwp_ax.pcolormesh(MEAStme_hgh_mesh,MEAShgt_hgh_mesh,MEASws_hgh,vmin = WSmin_hgh,\
        ws_rwp = rwp_ax.pcolormesh(MEAStme_hgh_mesh,MEAShgt_hgh_mesh,MEASws_hgh[:-1,:-1],vmin = WSmin_hgh,\
                                        vmax = WSmax_hgh,shading='flat',cmap=cmap_ws) 
        rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
        #######################################################################
        cbar = ws_fig.colorbar(ws_rwp,ax=rwp_ax,ticks = WStick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        # overlay barbs
        rwp_ax.barbs(MEAStme_hgh_mesh[BARBmeas_Hoi],MEAShgt_hgh_mesh[BARBmeas_Hoi],MEASu_hgh[BARBmeas_Hoi],\
                    MEASv_hgh[BARBmeas_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
        #######################################################################
        rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
        rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        # set figure y-label
        rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        # Denote the lowest measurement height
        rwp_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        #insert stn name/alt as title
        T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        ###############################

        #######################################################################
        # Plot the Model Wind Speed Time Histories
        ##########################################
        ws_mdl = mdl_ax.imshow(MEASws_hgh*np.nan,extent=RWPgrid_hgh,vmin=WSmin_hgh,vmax=WSmax_hgh,\
                             origin='lower',aspect='auto',cmap=cmap_ws);   
        mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        mdl_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); mdl_ax.set_yticks(HGTtick_hgh);mdl_ax.set_yticklabels(HGTtick_hgh_label);
        mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINhgh,MAXhgh]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        M = mdl_ax.set_title(MDLname + ' WS',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        M.set_position(plt_param.SPTLt_loc);
        
        """ ###################################################################
        Define and Plot the Model Wind Speed Biases
        ################################################################### """
        bias_ws = bias_ax.imshow(MEASws_hgh*np.nan,extent=RWPgrid_hgh,vmin = -3,vmax = 3,\
                             origin='lower',aspect='auto',cmap=cmap_bias);   
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
        #######################################################################
        bias_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); bias_ax.set_yticks(HGTtick_hgh);bias_ax.set_yticklabels(HGTtick_hgh_label);
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bias_ax.set_xlim(Xmin,Xmax)
        ####################################################################### 
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINhgh,MAXhgh]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        # Denote the lowest measurement height
        bias_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        #######################################################################
        B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
        ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WShgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        ws_fig.savefig(FILout,dpi=300); plt.close('all');
                             
        """ ###################################################################
        Wind Direction Plots (High [i.e. Low-Resolution] Mode)
        ################################################################### """
        wd_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                 height_ratios=[2,2,1],\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Lidar Wind Direction Time Histories
        ##############################################
        wd_fig.suptitle(str(frequency)+' MHz RWP WD',x = plt_param.DLmeas_x_title,\
                        y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
        #######################################################################
        wd_rwp = rwp_ax.pcolormesh(MEAStme_hgh_mesh,MEAShgt_hgh_mesh,MEASwd_hgh[:-1,:-1],vmin = WDmin_hgh,\
                                        vmax = WDmax_hgh,shading='flat',cmap=cmap_wd_hgh); 
        rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
        #######################################################################
        cbar = wd_fig.colorbar(wd_rwp,ax=rwp_ax,ticks = WDtick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        # overlay barbs
        rwp_ax.barbs(MEAStme_hgh_mesh[BARBmeas_Hoi],MEAShgt_hgh_mesh[BARBmeas_Hoi],MEASu_hgh[BARBmeas_Hoi],\
                    MEASv_hgh[BARBmeas_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
        #######################################################################
        rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
        rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        # set figure y-label
        rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        # Denote the lowest measurement height
        rwp_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        #insert stn name/alt as title
        T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        ###############################

        #######################################################################
        # Plot the Model Wind Direction Time Histories
        ##############################################
        wd_mdl = mdl_ax.imshow(MEASwd_hgh*np.nan,extent=RWPgrid_hgh,vmin=WDmin_hgh,vmax=WDmax_hgh,\
                             origin='lower',aspect='auto',cmap=cmap_wd_hgh); 
        mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        mdl_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); mdl_ax.set_yticks(HGTtick_hgh);mdl_ax.set_yticklabels(HGTtick_hgh_label);
        mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINhgh,MAXhgh]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        M = mdl_ax.set_title(MDLname + ' WD',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        M.set_position(plt_param.SPTLt_loc);
        
        """ ###################################################################
        Define and Plot the Model Wind Direction Biases
        ################################################################### """
        bias_wd = bias_ax.imshow(MEASwd_hgh*np.nan,extent=RWPgrid_hgh,vmin = -3,vmax = 3,\
                             origin='lower',aspect='auto',cmap=cmap_bias);   
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
        #######################################################################
        bias_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); bias_ax.set_yticks(HGTtick_hgh);bias_ax.set_yticklabels(HGTtick_hgh_label);
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINhgh,MAXhgh]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        # Denote the lowest measurement height
        bias_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
        wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WDhgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        wd_fig.savefig(FILout,dpi=300); plt.close('all');
                             
        """ ###################################################################
        Signal-to-Noise Ratio Plots (High [i.e. Low-Resolution] Mode)
        ################################################################### """
        #snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
        snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                 height_ratios=[3,1,1],\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Lidar Signal-to-Noise Time Histories
        ###############################################
        snr_fig.suptitle(str(frequency)+' MHz RWP Range-Corrected SNR',x = plt_param.DLmeas_x_title,\
                        y = plt_param.DLmeas_y_title,horizontalalignment='right',\
                fontsize=plt_param.TITLEfs);
        #######################################################################
        snr_rwp = rwp_ax.pcolormesh(MEAStme_hgh_mesh,MEAShgt_hgh_mesh,MEASsnr_hgh[:-1,:-1],vmin = SNRmin_hgh,\
                                        vmax = SNRmax_hgh,shading='flat',cmap=cmap_snr) 
        rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
        #######################################################################
        # overlay barbs
        rwp_ax.barbs(MEAStme_hgh_mesh[BARBmeas_Hoi],MEAShgt_hgh_mesh[BARBmeas_Hoi],MEASu_hgh[BARBmeas_Hoi],\
                    MEASv_hgh[BARBmeas_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
        #######################################################################
        cbar = snr_fig.colorbar(snr_rwp,ax=rwp_ax,ticks = SNRtick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect3)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('SNR, dB', fontsize=plt_param.LABELfs);
        #######################################################################
        rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
        rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        # set figure y-label
        rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        # Denote the lowest measurement height
        rwp_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        #insert stn name/alt as title
        T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        ###############################          
                             
        """ ###################################################################
        Plot the Model Wind Speed Biases
        ################################################################### """
        bias_ws = ws_ax.imshow(MEASws_hgh*np.nan,extent=RWPgrid_hgh,vmin = -3,vmax = 3,\
                             origin='lower',aspect='auto',cmap=cmap_ws);   
        ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(bias_ws,ax=ws_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
        #######################################################################
        ws_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); ws_ax.set_yticks(HGTtick_hgh);ws_ax.set_yticklabels(HGTtick_hgh_label);
        ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        ws_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        ws_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINhgh,MAXhgh]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        # Denote the lowest measurement height
        ws_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        B = ws_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################                     
                             
        """ ###################################################################
        Define and Plot the Model Wind Direction Biases
        ################################################################### """
        bias_wd = wd_ax.imshow(MEASwd_hgh*np.nan,extent=RWPgrid_hgh,vmin = -3,vmax = 3,\
                             origin='lower',aspect='auto',cmap=cmap_wd);   
        wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(bias_wd,ax=wd_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Dir, ($\mathrm{deg}^{\circ}$)', fontsize=plt_param.LABELfs);
        #######################################################################
        wd_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); wd_ax.set_yticks(HGTtick_hgh);wd_ax.set_yticklabels(HGTtick_hgh_label);
        wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        wd_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
        wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        #######################################################################
        wd_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([MINhgh,MAXhgh]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        wd_ax.set_xticks(Xtick); wd_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            wd_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        wd_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        wd_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        # Denote the lowest measurement height
        wd_ax.plot([Xmin,Xmax],[MEAShgt_hgh[0],MEAShgt_hgh[0]],'k--',linewidth=0.8)
        B = wd_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
        snr_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
        #######################################################################
        wd_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        wd_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_SNRhgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        snr_fig.savefig(FILout,dpi=300); plt.close('all');  
