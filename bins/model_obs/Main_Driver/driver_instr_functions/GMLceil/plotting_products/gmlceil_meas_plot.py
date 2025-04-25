#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Model Data Exists

###############################################################################
Created on Wed Jan 31 2023
@author: dmurray
###############################################################################
The Purpose of this Script is to Develop and Output time series plots
Plots Based on TROPoe retrievals
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Most Plots Generated in This Script Will Consist of Two Panels
- in the first plot, the obs will be plotted along with the model output
- in the second plot, the model bias will be shown and statistics will be annotated

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import os,time,warnings
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_mid import atmo_meas2hr_mid
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
##########################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
from avgmins_to_label import avgmins_to_label
from format_funcs import format_lat_lon
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
AVGper = getattr(plt_param,'CEIL_avgmins',30)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################

def gmlceil_meas(MEASin,OUTinfo):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth_model = os.path.join(OUTinfo[0],'no_model',DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    ###########################################################################
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=image.imread(plt_param.NOAAlogo)
    FILpfx='gmlceil_'
    ATTRstr = 'Credit: NOAA Global Monitoring Laboratory (GML)';
    
    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    gmlceil_plot(CEIL_plt,HRRRw_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""
    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat') 
    MEASalt = getattr(station_param,Sid+'lat')
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    
    ###########################################################################
    """ Unpack CEIL Information 
    ####################################
    MEASin :: CEIL_plt = (CEILloc,CEILtme,CEILhgt)
    ---> CEILloc = (CEILlon,CEILlat,CEILalt)
    ---> CEILtme = (CEILxy_minus,CEILxy,CEILxy_plus,CEILxy_ext)
    -----> CEILhgt = (CEIL_hgt_minus,CEIL_hgt,CEIL_hgt_plus,CEIL_hgt_ext)
    ###########################################################################
    CEILloc = (CEILlon,CEILlat,CEILalt); CEILxy = (time_meas,hgt);
    ###########################################################################
    CEIL_hgt = (hgt)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    CBHind_cat = -1; 
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEAShgt_in = MEASin[2]; 
    ###########################################################################
    MEAScbh = atmo_cat(MEAShgt_in[0],MEAShgt_in[1],MEAShgt_in[2],MEAShgt_in[3],CBHind_cat) 
    
    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    #MEAStime, MEAScbh = atmo_time_correct(MEAStime,MEAScbh);

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

        """ Define Relevant Plotting Specifics """
        # Define Proper Indices to Determine Axes Parameters
        meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
        ###############################################################
        # Perform Model Averaging for Plotting and Statistical Analsyes
        MEAScbh_mu = atmo_meas2hr_mid(MEAScbh,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
        ###############################################################
        CBHmin,CBHmax,CBHtick = atmo_bound_lim((MEAScbh[meas_oi]),\
                                       plt_param.CBHrnd_bse,plt_param.CBHmin,plt_param.CBHmax,\
                                       plt_param.CBHrnd_bse)
        """ ###############################################################
        Cloud Base Height Profile (instruct to share x-axis)
        ############################################################### """
        cbh_fig, (cbh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                 facecolor = 'w',edgecolor ='w',clear='True')
        ####################################################################
        # Plot the Cloud Base Height  Time Histories
        #####################################
        #Observations -- Only Plotting the QC'd Values
        cbh_ax.plot(MDLtime,MEAScbh_mu,'s',color='gray',markersize=plt_param.PLTmsize,label='Obs ('+AVGlabel+' Avg)')
        #cbh_ax.plot(MEAStime,MEAScbh, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
        cbh_ax.plot(MEAStime,MEAScbh,'.',color='black',markersize=1,linestyle='None',label='Obs')
        cbh_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        cbh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ####################################################################
        handles, labels = cbh_ax.get_legend_handles_labels()
        order = [1,0,2]
        cbh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        ###################################################################
        cbh_ax.set_ylim(CBHmin,CBHmax); cbh_ax.set_yticks(CBHtick);
        cbh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        ###################################################################
        cbh_ax.set_ylabel('Cloud Base Height,\nm AGL',{'fontsize': plt_param.LABELfs});
        cbh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        #insert station info as title
        T = cbh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);

        """ ###############################################################
        # Plot the Model Cloud Base Height Biases
        ############################################################### """
        # initate 0/0 line
        bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
        ###################################################################
        # Define the Bias (Observations - Measurements)
        bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
        bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ###################################################################
        bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                     ncol = 1,frameon = False);    
        ###################################################################
        bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
        bias_ax.set_xlim(Xmin,Xmax)
        ###################################################################
        bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
        bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);

        ####################################################################
        # Denote that Model Data was Not Available
        bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                     horizontalalignment='center',verticalalignment='bottom',color='red',\
                     fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        #Define Bounds of Precipitable Water - Model Bias to Determine Y-Bounds
        bias_ax.set_ylim(-100,100); bias_ax.set_yticks(np.arange(-100,100.1,50));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        bias_ax.set_ylabel('Model Error, m',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment': 'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        #######################################################################
        cbh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        ###################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        ###################################################################
        # add some attribution
        bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'center',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = cbh_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_CBH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        cbh_fig.savefig(FILout,dpi=300); plt.close('all');

