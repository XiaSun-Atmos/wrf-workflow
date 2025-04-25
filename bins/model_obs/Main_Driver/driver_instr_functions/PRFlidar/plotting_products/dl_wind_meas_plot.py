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
from datetime import datetime,timedelta,timezone
######################################
from time_tuple_cat import time_cat
from atmo_spatial_tuple_cat import atmo_spatial_cat
from dl_spatial_correction import dl_sptl_correct
###################################################
from atmo_ylim_plot import atmo_ylim
from atmo_spatial_ylim import atmo_sptl_ylim
from wd_spatial_ylim import wd_sptl_ylim
###################################################
from cmaps import get_var_cmap
from format_funcs import format_lat_lon_alt
####################################
#Initiate Plotting Parameter Module
import plotting_parameters as plt_param
import forecast_param_list as fxcst_param
import station_param_list as station_param

################################################################
# Define Model Interpolation Paramter
INTPoi = plt_param.INToi;
#0 indicates nearest-neighbor
#1 indicates bilinear interpolation
###############################################################################
HRsec = 60*60;
DAYsec = 24*HRsec;
MAXhgt = 200; #limit Doppler Lidar Plots to the Lowest 200 m

def dl_Wmeas(MEASin,OUTinfo):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
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
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    dl_Wplot(DLw_plt,HRRRw_plt,FILout)
    ###########################################################################
    MDLin ::  HRRRw_plt = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_wind)
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    ####################################################################### """
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
    """ Unpack Modeled Wind Data -- No Modelled Data Provided """
    #######################################################################

    #######################################################################
    # set up the colormaps
    cmap_ws = get_var_cmap('ws')
    cmap_bias = get_var_cmap('bias')
    cmap_wd = get_var_cmap('wd')

    ###########################################################################    
    # Define the Range of Possible Initialization Times #
    INIpos = ['%02d' % it for it in range(0,24)]
    """ Initiate Measurement Model Comparison for Each Initialization Hour """
    for ini in range(0,len(INIpos)):
        #######################################################################
        """ Initialization Time Not Available """
        #######################################################################
        MDLpseudo_time = np.arange(0,(FXlen-1)*3600+1e-5,step=3600)
        MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
        #######################################################################
        #Define Xtick Information (Only Needs to be Defined Once)
        Xmin = MDLtime[0];
        Xmax = (MDLtime[0]+(FXlen-1)*3600);
        #######################################################################
        Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
        Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
        Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
        Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
        Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
        
        #######################################################################
        # Define Proper Indices to Determine Axes Parameters
        meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[-1]))
        ###################################################################
        # pass through array indices only of interest here 
        WSmin,WSmax,WStick = atmo_sptl_ylim((MEASws[meas_oi]),\
                                plt_param.WSrnd_bse,plt_param.TICKint)    
        #WDmin,WDmax,WDtick = wd_sptl_ylim((MEASwd[meas_oi]));
        WDmin=0;WDmax=360;WDtick=np.arange(0,361,60);
        ###################################################################
        # maximum and minimum measurement height 
        HGTmin,HGTmax,HGTtick = atmo_ylim((np.asarray([0,MAXhgt])),plt_param.SHGTrnd_bse,plt_param.TICKint);
        HGTtick_label = [str(int(ht)) for ht in HGTtick]
        
        #######################################################################
        # establish barb measurement array
        BARBmeas_oi = np.zeros(MEAStme_mesh.shape,dtype = bool)
        #######################################################################
        # define the time stamps of interest
        BARBtme_oi = (60*60)*np.arange(ini + plt_param.BARBdl_x_base,(ini + FXlen),\
                      plt_param.BARBdl_x);
        #######################################################################
        # Define the Max Height Ind
        MEAShgt_max = np.where(MEAShgt > MAXhgt)[0][0] 
        #######################################################################
        for j in np.arange(plt_param.BARBdl_y_base,MEAShgt_max,plt_param.BARBdl_y):
            for toi in BARBtme_oi:
                # identify the forecast time reference indice #
                barb_meas = np.where(np.abs(MEAStime - toi) - np.min(np.abs(MEAStime - toi)) \
                                     < 1e-5)[0][0]
                ###############################################################
                # only keep a single indice if equi-distant
                ###############################################################
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
        ws_mdl = mdl_ax.imshow(MEASws*np.nan,extent=DLgrid,vmin=WSmin,vmax=WSmax,\
                             origin='lower',aspect='auto', cmap = cmap_ws);   
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
        # plotting function will auto-expand the grid
        ###################################################################
        M = mdl_ax.set_title(MDLname + ' WS',\
                             {'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        M.set_position(plt_param.SPTLt_loc);

        """ ###############################################################
        Define and Plot the Model Wind Speed Biases
        ############################################################### """
        # Plot the Model Wind Speed Bias
        #######################################################################
        bias_ws = bias_ax.imshow(MEASws*np.nan,extent=DLgrid,vmin=-3,vmax=3,\
                             origin='lower',aspect='auto', cmap = cmap_bias);   
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
        #######################################################################
        bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        bias_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        ###################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        # Denote that Model Data was Not Available
        bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
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
        bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'center',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
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
        ##########################################
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
        ###################################################################

        ###################################################################
        # Plot the Model Wind Direction Time Histories
        ##########################################
        wd_mdl = mdl_ax.imshow(MEASwd*np.nan,extent=DLgrid,vmin=WDmin,vmax=WDmax,\
                             origin='lower',aspect='auto', cmap = cmap_wd);   
        mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ###################################################################
        cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
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
        # plotting function will auto-expand the grid
        ###################################################################
        M = mdl_ax.set_title(MDLname + ' WD',\
                             {'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        M.set_position(plt_param.SPTLt_loc);

        """ ###############################################################
        Define and Plot the Model Wind Direction Biases
        ############################################################### """
        # Plot the Model Wind Direction Bias
        #######################################################################
        bias_wd = bias_ax.imshow(MEASwd*np.nan,extent=DLgrid,vmin=-3,vmax=3,\
                             origin='lower',aspect='auto', cmap = cmap_bias);   
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
        cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
        cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
        #######################################################################
        bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
        #######################################################################
        # set figure y-label
        bias_ax.set_ylabel('Height, m AGL',{'fontsize': plt_param.LABELfs});
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
        ###################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=UserWarning)
            bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        # Denote that Model Data was Not Available
        bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                 fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                     horizontalalignment='right');  
        B.set_position(plt_param.SPTLt_loc);
        #######################################################################
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
        bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'center',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = wd_fig.add_axes(plt_param.LOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################          
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        wd_fig.savefig(FILout,dpi=300); plt.close('all');
                             
