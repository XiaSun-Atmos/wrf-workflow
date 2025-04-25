#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case where both Measurement and Model Data are 
Available

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Temp/Pressure/Mixing Ratio
Plots Based on Measurements from the SFC Standard Meteorological Information 
Installed at Various Sites
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Most Plots Generated in This Script Will Consist of Two Panels 
- in the first plot, the obs will be plotted along with the model output
- in the second plot, the model bias will be shown and statistics will be annotated

- Model Data Still Needs to be Incorporated

This Function Will No Longer Output Plots of Relative Humidity -- Instead Plots
of the Mixing Ratio will be Produced
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import os,time,warnings
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_end import atmo_meas2hr_end
from atmo_meas_to_hour_upto import atmo_meas2hr_upto
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
##########################################
from atmo_bias import atmo_bias
from atmo_ylim_bias import atmo_blim
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
from label_funcs import avgmins_to_label
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
# averaging period for comparison (minutes), label for plot
AVGper = getattr(plt_param,'PSLmet_avgmins',60)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################

def psl_RADmeas(MEASin,OUTinfo):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth = os.path.join(OUTinfo[0],DATE);
    if not os.path.exists(IMGpth):
        os.makedirs(IMGpth)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext',DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_model = os.path.join(OUTinfo[0],'no_model',DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    ###########################################################################Z
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=image.imread(plt_param.NOAAlogo);
    FILpfx='pslmet_'
    
    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #    EXTlen = fxcst_param.hrrr_v4_ext;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    #    EXTlen = fxcst_param.rap_v5_ext;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    psl_RADmeas(MEAStrp_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""

    ###########################################################################
    """ No Model Data to Unpack
    ###############################################################"""
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
    """ Unpack Measured RAD Information 
    ####################################
    MEASrad_plt = (MEASloc,MEAStme,MEASrad_in)       
    ---> MEASloc = (MEASlon,MEASlat,MEASalt)
    ---> MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus)
    ---> MEASrad_in = (MEASrad_minus,MEASrad,MEAStrh_plus)
    ###########################################################################
    MEASrad = (MEASsolar,MEASnet_rad)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    SRADind_cat = 0; NRADind_cat = 1; 
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASrad_in = MEASin[2]
    ######################
    MEASsolar = atmo_cat(MEASrad_in[0],MEASrad_in[1],MEASrad_in[2],MEASrad_in[3],SRADind_cat) 
    MEASsolar_qc = None

    MEASnet_rad = atmo_cat(MEASrad_in[0],MEASrad_in[1],MEASrad_in[2],MEASrad_in[3],NRADind_cat) 
    MEASnet_rad_qc = None

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStime_new, MEASsolar = atmo_time_correct(MEAStime,MEASsolar);
    _, MEASnet_rad = atmo_time_correct(MEAStime,MEASnet_rad);
    MEAStime = MEAStime_new

    # Define the Range of Possible Initialization Times #
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

        ###############################################################
        # Perform Model Averaging for Plotting and Statistical Analsyes
        MEASsolar_mu = atmo_meas2hr_end(MEASsolar,MEASsolar_qc,MEAStime,MDLtime, AVGmins=AVGper);
        MEASnet_rad_mu = atmo_meas2hr_end(MEASnet_rad,MEASnet_rad_qc,MEAStime,MDLtime, AVGmins=AVGper);

        """ Define Relevant Plotting Specifics """
        # Define Proper Indices to Determine Axes Parameters
        meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
        ###############################################################
        SRADmin,SRADmax,SRADtick = atmo_ylim((MEASsolar[meas_oi]),\
                                       plt_param.SRADrnd_bse,plt_param.TICKint)
        NRADmin,NRADmax,NRADtick = atmo_ylim((MEASnet_rad[meas_oi]),\
                                       plt_param.SRADrnd_bse,plt_param.TICKint)

        ###################################################################
        """ Measurements and Model Data are Available """
        """ ###############################################################
        Shortwave Radiation Measurements
        ############################################################### """
        rad_fig, (sw_ax,net_ax,bias_ax) = plt.subplots(3,1,sharex = True, \
                  figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),facecolor = 'w',edgecolor ='w',clear='True')
        ###################################################################
        # Plot the Shortwave Downwelling Time Histories
        ###############################################
        sw_ax.plot(MDLtime,MEASsolar_mu,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs (1Hr Avg)')
        sw_ax.plot(MEAStime,MEASsolar,'-',color='black',linewidth=plt_param.PLTlwidth,markersize=plt_param.PLTmsize,label='Obs')
        sw_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        sw_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ####################################################################
        handles, labels = sw_ax.get_legend_handles_labels()
        order = [1,0,2]
        sw_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       #loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_tri_ur,\
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = (1.025,1.25),\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        ###################################################################
        sw_ax.set_ylim(SRADmin,SRADmax); sw_ax.set_yticks(SRADtick);
        sw_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        ###################################################################
        sw_ax.set_ylabel('Solar Radiation,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});
        sw_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        #insert station info as title
        T = sw_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
                     
        ###################################################################
        # Plot the Net Radiation Time Histories
        ###############################################
        net_ax.plot(MDLtime,MEASnet_rad_mu,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs (1Hr Avg)')
        net_ax.plot(MEAStime,MEASnet_rad,'-',color='black',linewidth=plt_param.PLTlwidth,markersize=plt_param.PLTmsize,label='Obs')
        net_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        net_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ####################################################################
        net_ax.set_ylim(NRADmin,NRADmax); net_ax.set_yticks(NRADtick);
        net_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        ###################################################################
        net_ax.set_ylabel('Net Radiation,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});
        net_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        # 0/0 line
        net_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=.5)
    
        """ ###############################################################
        # Plot the Model Radiation Biases
        ############################################################### """
        # initate 0/0 line
        bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
        bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='forestgreen',markersize=plt_param.PLTmsize,\
             linewidth=plt_param.PLTlwidth,label='Solar Radiation (RMSE: N/A)')
        bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='dodgerblue',markersize=plt_param.PLTmsize,\
             linewidth=plt_param.PLTlwidth,label='Net Radiation (RMSE: N/A)')
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ###################################################################
        #bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
        #             ncol = 1,frameon = False);  
        #######################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        bias_ax.set_xlim(Xmin,Xmax)
        #######################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
        #######################################################################
        # Denote that Model Data was Not Availabel
        bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                     horizontalalignment='center',verticalalignment='bottom',color='red',\
                     fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
        bias_ax.set_ylim(-100,100); bias_ax.set_yticks(np.arange(-100,100.1,50));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        ####################################################################
        bias_ax.set_ylabel('Model Error,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        bias_ax.legend(loc = plt_param.BIASpad,fontsize=plt_param.LEGfs,ncol = 2,\
                              frameon = False);
        ###################################################################
        bias_ax.annotate('Mdl - Obs',xy=(0.13,0.345),xycoords='figure fraction',color='black',\
                     fontsize=plt_param.LEGfs,weight='bold');
        ###################################################################
        rad_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        ###################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = rad_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth,FILpfx+MDLname+'_SRAD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        rad_fig.savefig(FILout,dpi=300); plt.close('all');
    
