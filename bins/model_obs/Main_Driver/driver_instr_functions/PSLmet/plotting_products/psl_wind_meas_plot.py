#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Model Data Exists
    
###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
This Script is a Variant of psl_wind_plot. However, in this script, only 
measurement data were defined so that is the only variable to be plotted
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Most Plots Generated in This Script Will Consist of Two Panels 
- in the first plot, the obs will be plotted along with the model output
- in the second plot, the model bias will be shown and statistics will be annotated

###############################################################################
    #For Statisical Wind Direction Computations -- Only Consider Wind Speeds 
    #Abovce Some Threshold Value

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import os,time
from datetime import datetime,timedelta
#####################################
from atmo_meas_to_hour_end import atmo_meas2hr_end
from wd_meas_to_hour_end import wd_meas2hr_end
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
#####################################
from atmo_ylim_plot import atmo_ylim
from wd_ylim_plot import wd_ylim
from wd_disc import wd_plt
from avgmins_to_label import avgmins_to_label
from format_funcs import format_lat_lon
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
AVGper = getattr(plt_param,'PSLmet_avgmins',60)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################
# Define SFC Measurement Height 
SFCws_hgt = 10;
# per Dan Gottas -- For the met systems: winds = 10m, pressure = 2m, temp/rh = 2m.


def psl_Wmeas(MEASin,OUTinfo):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth_model = os.path.join(OUTinfo[0],'no_model',DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    ###########################################################################
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=img.imread(plt_param.NOAAlogo)
    
    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    psl_Wmeas(MEASin,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""

    ###########################################################################
    """ Modeled Wind Data Not Available - No Unpacking Necessary """
    ###########################################################################
        
    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat') 
    MEASalt = getattr(station_param,Sid+'alt')
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    ###########################################################################
    
    ###########################################################################
    """ Unpack Measured Wind Information 
    ####################################
    MEASw_plt = (MEASloc,MEAStme,MEASw_in)         
    ---> MEASloc = (MEASlon,MEASlat,MEASalt)
    ---> MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus)
    ---> MEASw_in = (MEASwind_minus,MEASwind,MEASwind_plus)
    ###########################################################################
    MEASw = (MEASws,MEASwd_vec,MEASwd_vec_std)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    WSind_cat = 0; WSqc_ind_cat = 1;
    WDind_cat = 1; WDqc_ind_cat = 4;
    ###########################################################################
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASw_in = MEASin[2]; 
    #####################
    MEASws = atmo_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WSind_cat)
    MEASws_qc = None
    MEASwd = atmo_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WDind_cat)
    MEASwd_qc = None

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStime_new, MEASws = atmo_time_correct(MEAStime,MEASws);
    _, MEASwd = atmo_time_correct(MEAStime,MEASwd);
    MEAStime = MEAStime_new;

    # Define the Range of Possible Initialization Times to be Considered #
    INIpos = ['%02d' % it for it in range(0,24)]
    """ Initiate Measurement Model Comparison for Each Initialization Hour """
    for ini in range(0,len(INIpos)):
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
                
        """ Perform Model Averaging for Plotting and Statistical Analyses """
        MEASws_mu = atmo_meas2hr_end(MEASws,MEASws_qc,MEAStime,MDLtime, AVGmins=AVGper);
        MEASwd_mu = wd_meas2hr_end(MEASwd,MEASwd_qc,MEAStime,MDLtime, AVGmins=AVGper);
        ####################################################################### 
        
        """ Define Relevant Plotting Specifics """
        # Define Proper Indices to Determine Axes Parameters
        meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
        #######################################################################
        #WSmin,WSmax,WStick = atmo_ylim((MEASws[np.logical_and(MEASws_qc == 0,meas_oi)]),\
        WSmin,WSmax,WStick = atmo_ylim((MEASws[meas_oi]),\
                                       plt_param.WSrnd_bse,plt_param.TICKint)
        #WDmin,WDmax,WDtick = wd_ylim((MEASwd[np.logical_and(MEASwd_qc == 0,meas_oi)]));
        WDmin,WDmax,WDtick = wd_ylim((MEASwd[meas_oi]));

        """ ###################################################################
        Wind Speed Plots (instruct to share x-axis)
        ################################################################### """
        ws_fig, (ws_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Wind Speed Time Histories
        ####################################
        #Observations -- Only Plotting the QC'd Values
        ws_ax.plot(MDLtime,MEASws_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
        ws_ax.plot(MEAStime,MEASws,color='black', linewidth=plt_param.PLTlwidth,label='Obs')
        ws_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        handles, labels = ws_ax.get_legend_handles_labels()
        order = [1,0,2]
        ws_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        ws_ax.set_ylim(WSmin,WSmax); ws_ax.set_yticks(WStick);
        ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        ws_ax.set_ylabel('Wind Speed, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});
        ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert location as title
        T = ws_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
                
        """ ###################################################################
        Define and Plot the Model Wind Speed Biases
        #################################################################### """
        # initate 0/0 line
        bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
        #######################################################################
        bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ###################################################################
        bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                     ncol = 1,frameon = False);  
        #######################################################################
        #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
        #######################################################################
        # Denote that Model Data was Not Availabel
        bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                     horizontalalignment='center',verticalalignment='bottom',color='red',\
                     fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
        bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        bias_ax.set_ylabel('Model Error, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        #######################################################################
        #-> annotate measurement heights
        ws_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                         fontsize=plt_param.LEGfs,weight='bold')
        #######################################################################
        ws_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,'pslmet_'+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        ws_fig.savefig(FILout,dpi=300); plt.close('all');
                
        """ ###################################################################
        Wind Direction Plots (instruct to share x-axis)
        ################################################################### """
        wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                 facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Wind Direction Time Histories
        #######################################################################
        #MEASwd[MEASwd_qc != 0] = np.nan;
        wd_ax = wd_plt(MDLtime,MEASwd_mu,wd_ax,'s','gray','Obs (1Hr Avg)',marker_size=plt_param.PLTmsize)
        wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs')
        wd_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        #######################################################################
        wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        wd_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                     ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);             
        #######################################################################
        wd_ax.set_ylim(WDmin,WDmax); wd_ax.set_yticks(WDtick);
        wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        wd_ax.set_ylabel('Wind Direction, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});
        wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###########################################################################
        #insert location as title
        T = wd_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
                
        """ ###################################################################
        # Plot the Model Wind Direction Biases
        ################################################################### """
        # initate 0/0 line
        bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
        #######################################################################
        bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)   
        ###################################################################
        bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                     ncol = 1,frameon = False);  
        #######################################################################
        #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
        #######################################################################
        # Denote that Model Data was Not Availabel
        bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                     horizontalalignment='center',verticalalignment='bottom',color='red',\
                     fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
        bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        bias_ax.set_ylabel('Model Error, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        #######################################################################
        #-> annotate measurement heights
        wd_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                         fontsize=plt_param.LEGfs,weight='bold')
        #######################################################################
        wd_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
                        #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,'pslmet_'+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        wd_fig.savefig(FILout,dpi=300); plt.close('all')
        #######################################################################
