#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Model Data Exists

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Precipitation
Plots Based on Measurements from the PSL Standard Meteorological Information 
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
###############################################################################
This Function Will No Longer Output Plots of Relative Humidity -- Instead Plots
of the Mixing Ratio will be Produced
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import pandas as pd
from datetime import datetime,timedelta
import os,time,warnings
#####################################
from atmo_meas_to_hour_end import atmo_meas2hr_end
from atmo_meas_to_hour_beg import atmo_meas2hr_beg
from atmo_meas_to_hour_upto import atmo_meas2hr_upto
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
##########################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
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
###############################################################################

def psl_PRmeas(MEASin,OUTinfo):
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
    FILpfx='psldisd_'
    
    """ Based on MDLname -- Define Relevant Variables """
    if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
        FXlen = fxcst_param.hrrr_v4_fx;
    elif MDLname == 'RAP_v5': # 'ESRL_RAP';
        FXlen = fxcst_param.rap_v5_fx;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    psl_TRPplot(MEAStrp_plt,HRRRtrp_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""
    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat')
    MEASalt = getattr(station_param,Sid+'alt')
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    
    ###########################################################################
    """ Unpack Measured TRP Information 
    ####################################
    MEAStrp_plt = (MEASloc,MEAStme,MEAStrh_in,MEASpres_in)       
    ---> MEASloc = (MEASlon,MEASlat,MEASalt)
    ---> MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus)
    ---> MEAStrh_in = (MEAStrh_minus,MEAStrh,MEAStrh_plus)
    ---> MEASpres_in = (MEASpres_minus,MEASpres,MEASpres_plus)
    ###########################################################################
    MEAStrh = (MEAStemp,MEASrh)
    MEASpres = (MEASpress)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    PRATEind_cat = 0; 
    PRind_cat = 1; 
    PRSUMind_cat = 2;
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASpr_in = MEASin[2]
    ######################
    MEASpr = atmo_cat(MEASpr_in[0],MEASpr_in[1],MEASpr_in[2],MEASpr_in[3],PRind_cat) 

    # Define the Range of Possible Initialization Times #
    INIpos = ['%02d' % it for it in range(0,24)]
    """ Initiate Measurement Model Comparison for Each Initialization Hour """
    for ini in range(0,len(INIpos)):
        #######################################################################
        MDLpseudo_time = np.arange(0,(FXlen-1)*3600+1e-5,step=3600)
        MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
        #######################################################################
        #Define Xtick Information (Only Needs to be Defined Once)
        Xmin = MDLtime[0] - 3600;
        Xmax = (MDLtime[0]+(FXlen-1)*3600) + 3600;
        #######################################################################
        Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
        Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
        Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
        Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
        Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')

        """ Perform Model Averaging for Plotting and Statistical Analyses """
        MEASpr_mu = atmo_meas2hr_upto(MEASpr,None,MEAStime,MDLtime,func='sum');
        #######################################################################
        
        """ Define Relevant Plotting Specifics """
        # Define Proper Indices to Determine Axes Parameters
        meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
        #######################################################################
        PRmin,PRmax,PRtick = atmo_bound_lim((MEASpr[meas_oi]),\
                     plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                     plt_param.PRrnd_bse)
        APCPmin,APCPmax,APCPtick = atmo_bound_lim((MEASpr_mu),\
                     plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                     plt_param.APCPrnd_bse)

        """ ###################################################################
        Precipitation Profile (instruct to share x-axis)
        ################################################################### """
        pr_fig, (apcp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                   facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Precipitation Time Histories
        #####################################
        #Observations -- Only Plotting the QC'd Values
        pr_ax = apcp_ax.twinx()
        pr_ax.bar(MEAStime,MEASpr,width=120,color='blue',label='Obs',align='edge')
        
        apcp_ax.plot(MDLtime,MEASpr_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='1-Hour Total')
        apcp_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        apcp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        prlines,prlabels = pr_ax.get_legend_handles_labels()
        apcplines,apcplabels = apcp_ax.get_legend_handles_labels()
        apcp_ax.legend(prlines+apcplines, prlabels+apcplabels, loc=1,\
                        fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                        ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        apcp_ax.set_ylim(APCPmin,APCPmax); apcp_ax.set_yticks(APCPtick);
        apcp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        pr_ax.set_ylim(PRmin,PRmax); pr_ax.set_yticks(PRtick);
        pr_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        apcp_ax.set_ylabel('1-Hr Total Precip, mm',{'fontsize': plt_param.LABELfs});
        apcp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        pr_ax.set_ylabel('2-min Obs Precip, mm',{'fontsize': plt_param.LABELfs});
        #######################################################################
        #insert station info as title
        T = apcp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        
        """ ###################################################################
        Plot the Model Precipitation Biases
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
        bias_ax.set_ylim(-5,5); bias_ax.set_yticks(np.arange(-5,5.1,2.5));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        ###################################################################
        #######################################################################
        pr_fig.subplots_adjust(left = 0.125,right = 0.905,bottom = 0.155, top = 0.935)
        #######################################################################
        bias_ax.annotate(DATE[0:4]+'-'+DATE[4:6]+'-'+DATE[6:9],xy = plt_param.STNpad_sptl,xycoords = 'figure fraction',\
                verticalalignment = 'bottom',horizontalalignment = 'center',\
                fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = pr_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_PR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        pr_fig.savefig(FILout,dpi=300); plt.close('all');
        #######################################################################
        
