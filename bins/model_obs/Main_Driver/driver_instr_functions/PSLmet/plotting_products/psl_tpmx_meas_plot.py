#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Model Data Exists

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Temp/Pressure/Mixing Ratio
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
import os,time,warnings
from datetime import datetime,timedelta
#####################################
from atmo_meas_to_hour_end import atmo_meas2hr_end
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
from atmo_moisture import calc_dewpoint_from_RH
from atmo_moisture import calc_mixingratio
##########################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
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
# averaging period for comparison (minutes), label for plot
AVGper = getattr(plt_param,'PSLmet_avgmins',60)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################
# Define PSL Measurement Heights 
SFCpress_hgt = 1;
SFCtemp_hgt = 2;
SFCrh_hgt = 2;
# per Jenni Kyrouac -- For the met systems: winds = 10m, pressure = 1m, temp/rh = 2m.

def psl_TPMXmeas(MEASin,OUTinfo):
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
    FILpfx='pslmet_'
    
    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
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
    MEASalt = getattr(station_param,Sid+'lat')
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
    ---> MEASprecip_in = (MEASprecip_minus,MEASprecip,MEASprecip_plus)
    ###########################################################################
    MEAStrh = (MEAStemp,MEASrh)
    MEASpres = (MEASpress)
    MEASpr = (MEASprecip)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    TEMPind_cat = 0; RHind_cat = 1; 
    PRESind_cat = -1; 
    PRind_cat = -1; 
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEAStrh_in = MEASin[2]
    MEASprs_in = MEASin[3]
    MEASpr_in = MEASin[4]
    ######################
    MEAStemp = atmo_cat(MEAStrh_in[0],MEAStrh_in[1],MEAStrh_in[2],MEAStrh_in[3],TEMPind_cat) 
    MEAStemp_qc = None
    
    MEASrh = atmo_cat(MEAStrh_in[0],MEAStrh_in[1],MEAStrh_in[2],MEAStrh_in[3],RHind_cat) 
    MEASrh_qc = None
    
    MEASpres = atmo_cat(MEASprs_in[0],MEASprs_in[1],MEASprs_in[2],MEASprs_in[3],PRESind_cat); 
    MEASpres_qc = None
    
    MEASprecip = atmo_cat(MEASpr_in[0],MEASpr_in[1],MEASpr_in[2],MEASpr_in[3],PRind_cat) 

    MEAStd = calc_dewpoint_from_RH(MEAStemp,MEASrh)
    MEASlcl = 125*(MEAStemp-MEAStd)
    
    MEASmxrat = calc_mixingratio(MEAStemp,MEASrh, MEASpres)
    MEASmxrat_qc = None

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStime_new, MEAStemp = atmo_time_correct(MEAStime,MEAStemp);
    _, MEASrh = atmo_time_correct(MEAStime,MEASrh);
    _, MEASpres = atmo_time_correct(MEAStime,MEASpres);
    _, MEASprecip = atmo_time_correct(MEAStime,MEASprecip);
    _, MEASlcl = atmo_time_correct(MEAStime,MEASlcl);
    _, MEASmxrat = atmo_time_correct(MEAStime,MEASmxrat);
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

        """ Perform Model Averaging for Plotting and Statistical Analyses """
        MEAStemp_mu = atmo_meas2hr_end(MEAStemp,MEAStemp_qc,MEAStime,MDLtime, AVGmins=AVGper);
        MEASpres_mu = atmo_meas2hr_end(MEASpres,MEASpres_qc,MEAStime,MDLtime, AVGmins=AVGper);
        MEASmxrat_mu = atmo_meas2hr_end(MEASmxrat,MEASmxrat_qc,MEAStime,MDLtime, AVGmins=AVGper)
        MEASprecip_mu = atmo_meas2hr_end(MEASprecip,None,MEAStime,MDLtime, AVGmins=AVGper)
        MEASrh_mu = atmo_meas2hr_end(MEASrh,MEASmxrat_qc,MEAStime,MDLtime, AVGmins=AVGper)
        MEASlcl_mu = atmo_meas2hr_end(MEASlcl,None,MEAStime,MDLtime, AVGmins=AVGper)
        #######################################################################
        
        """ Define Relevant Plotting Specifics """
        # Define Proper Indices to Determine Axes Parameters
        meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
        #######################################################################
        TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MEAStemp[meas_oi]),\
                                       plt_param.TEMPrnd_bse,plt_param.TICKint)
        PRESmin,PRESmax,PREStick = atmo_ylim((MEASpres[meas_oi]),\
                                       plt_param.PRESrnd_bse,plt_param.TICKint)
        MXRATmin,MXRATmax,MXRATtick = atmo_ylim((MEASmxrat[meas_oi]),\
                                       plt_param.MXRATrnd_bse,plt_param.TICKint)
        #RHmin,RHmax,RHtick = atmo_ylim((MEASrh[meas_oi]),\
        #                               plt_param.RHrnd_bse,plt_param.TICKint)
        RHmin = 0; RHmax = 110; RHtick = np.arange(RHmin,RHmax+1,20)
        LCLmin,LCLmax,LCLtick = atmo_ylim((MEASlcl[meas_oi]),\
                                       plt_param.LCLrnd_bse,plt_param.TICKint)
        LCLmin = 0;LCLmax = max(1500,LCLmax); LCLtick = np.arange(LCLmin,LCLmax+1,plt_param.LCLrnd_bse);
        PRmin,PRmax,PRtick = atmo_bound_lim((MEASprecip[meas_oi]),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
        APCPmin,APCPmax,APCPtick = atmo_bound_lim((MEASprecip_mu),\
                                 plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
        """ ###################################################################
        Temperature Profile (instruct to share x-axis)
        ################################################################### """
        temp_fig, (temp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                   facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Temperature Time Histories
        #####################################
        #Observations 
        temp_ax.plot(MDLtime,MEAStemp_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
        temp_ax.plot(MEAStime,MEAStemp, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
        temp_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        handles, labels = temp_ax.get_legend_handles_labels()
        order = [1,0,2]
        temp_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        temp_ax.set_ylim(TEMPmin,TEMPmax); temp_ax.set_yticks(TEMPtick);
        temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        temp_ax.set_ylabel('Temperature, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});
        temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert station info as title
        T = temp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        
        """ ###################################################################
        Plot the Model Temperature Biases
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
        bias_ax.set_ylabel('Model Error, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        ###################################################################
        #-> annotate measurement heights
        temp_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                         fontsize=plt_param.LEGfs,weight='bold')
        #######################################################################
        temp_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = temp_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        temp_fig.savefig(FILout,dpi=300); plt.close('all');
        #######################################################################
        
        """ ###################################################################
        Pressure Plots (instruct to share x-axis)
        ################################################################### """
        pres_fig, (pres_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                   facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Atmospheric Pressure Time Histories
        ##############################################
        #Observations 
        pres_ax.plot(MDLtime,MEASpres_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
        pres_ax.plot(MEAStime,MEASpres,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
        pres_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        pres_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        handles, labels = pres_ax.get_legend_handles_labels()
        order = [1,0,2]
        pres_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        pres_ax.set_ylim(PRESmin,PRESmax); pres_ax.set_yticks(PREStick);
        pres_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        pres_ax.set_ylabel('Pressure, mb',{'fontsize': plt_param.LABELfs});
        pres_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert station info as title
        T = pres_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        
        """ ###################################################################
        # Plot the Model Pressure Biases
        ##################################################################  """
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
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        #######################################################################
        bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        #-> annotate measurement heights
        pres_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                         fontsize=plt_param.LEGfs,weight='bold')
        #######################################################################
        pres_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = pres_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_PRES_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        pres_fig.savefig(FILout,dpi=300); plt.close('all');

        """ ###################################################################
        Lifting Condenstation Level Plots (instruct to share x-axis)
        ################################################################### """
        lcl_fig, (lcl_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                   facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Atmospheric Pressure Time Histories
        ##############################################
        #Observations 
        lcl_ax.plot(MDLtime,MEASlcl_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
        lcl_ax.plot(MEAStime,MEASlcl,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
        lcl_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        lcl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        handles, labels = lcl_ax.get_legend_handles_labels()
        order = [1,0,2]
        lcl_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        lcl_ax.set_ylim(LCLmin,LCLmax); lcl_ax.set_yticks(LCLtick);
        lcl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        lcl_ax.set_ylabel('Lifting Condensation \nLevel, m',{'fontsize': plt_param.LABELfs});
        lcl_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert station info as title
        T = lcl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        
        """ ###################################################################
        # Plot the Model Pressure Biases
        ##################################################################  """
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
        bias_ax.set_ylim(-250,250); bias_ax.set_yticks(np.arange(-250,251,250));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        #######################################################################
        bias_ax.set_ylabel('Model Error, m',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        lcl_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = lcl_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_LCL_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        lcl_fig.savefig(FILout,dpi=300); plt.close('all');

        """ ###################################################################
        Mixing Ratio Profile (instruct to share x-axis)
        ################################################################### """
        mxrat_fig, (mxrat_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                   facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Mixing Ratio Time Histories
        #####################################
        #Observations 
        mxrat_ax.plot(MDLtime,MEASmxrat_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
        mxrat_ax.plot(MEAStime,MEASmxrat, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
        mxrat_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        mxrat_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        handles, labels = mxrat_ax.get_legend_handles_labels()
        order = [1,0,2]
        mxrat_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        mxrat_ax.set_ylim(MXRATmin,MXRATmax); mxrat_ax.set_yticks(MXRATtick);
        mxrat_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        mxrat_ax.set_ylabel('Mixing Ratio, g kg-1',{'fontsize': plt_param.LABELfs});
        mxrat_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert station info as title
        T = mxrat_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        
        """ ###################################################################
        Plot the Model Mixing Ratio Biases
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
        bias_ax.set_ylabel('Model Error, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        ###################################################################
        #-> annotate measurement heights
        mxrat_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                         fontsize=plt_param.LEGfs,weight='bold')
        #######################################################################
        mxrat_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = mxrat_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_MXRAT_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        mxrat_fig.savefig(FILout,dpi=300); plt.close('all');
        #######################################################################
        
        """ ###################################################################
        Relative Humidity Plots (instruct to share x-axis)
        ################################################################### """
        rh_fig, (rh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                   facecolor = 'w',edgecolor ='w',clear='True')
        #######################################################################
        # Plot the Atmospheric Relative Humidity Time Histories
        ##############################################
        #Observations 
        rh_ax.plot(MDLtime,MEASrh_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
        rh_ax.plot(MEAStime,MEASrh,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
        rh_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        rh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        #######################################################################
        handles, labels = rh_ax.get_legend_handles_labels()
        order = [1,0,2]
        rh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                       loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        #######################################################################
        rh_ax.set_ylim(RHmin,RHmax); rh_ax.set_yticks(RHtick);
        rh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        rh_ax.set_ylabel('Relative Humidity, %',{'fontsize': plt_param.LABELfs});
        rh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        #######################################################################
        #insert station info as title
        T = rh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
        
        """ ###################################################################
        # Plot the Model Relative Humidity Biases
        ##################################################################  """
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
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        #######################################################################
        bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        #-> annotate measurement heights
        rh_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                         fontsize=plt_param.LEGfs,weight='bold')
        #######################################################################
        rh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
        #######################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        #######################################################################
        # add a logo
        imax = rh_fig.add_axes(plt_param.metLOGOaxcoords)
        imax.set_axis_off()
        imax.imshow(logo, aspect="equal")
        #######################################################################
        #Save the Developed Figure
        FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_RH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
        rh_fig.savefig(FILout,dpi=300); plt.close('all');

        """ ###############################################################
        Precipitation Profile (instruct to share x-axis)
        ############################################################### """
        pr_fig, (apcp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                 facecolor = 'w',edgecolor ='w',clear='True')
        ####################################################################
        # Plot the Precipitation Time Histories
        #####################################
        #Observations -- Only Plotting the QC'd Values
        pr_ax=apcp_ax.twinx();
        pr_ax.bar(MEAStime,MEASprecip,width=120,color='blue',label='Obs',align='edge',linewidth=0.5)
        
        apcp_ax.plot(MDLtime,MEASprecip_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='1-Hour Total')
        apcp_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
        apcp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
        ####################################################################
        prlines,prlabels = pr_ax.get_legend_handles_labels()
        apcplines,apcplabels = apcp_ax.get_legend_handles_labels()
        apcp_ax.legend(prlines+apcplines, prlabels+apcplabels, loc=1,\
                       fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                       ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
        ###################################################################
        apcp_ax.set_ylim(APCPmin,APCPmax); apcp_ax.set_yticks(APCPtick);
        apcp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        pr_ax.set_ylim(PRmin,PRmax); pr_ax.set_yticks(PRtick);
        pr_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        ###################################################################
        apcp_ax.set_ylabel('1-Hr Total Precip, mm',{'fontsize': plt_param.LABELfs});
        apcp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        pr_ax.set_ylabel('2-min Obs Precip, mm',{'fontsize': plt_param.LABELfs});
        ###################################################################
        #insert station info as title
        T = apcp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        T.set_position(plt_param.Tloc);
                     
        """ ###############################################################
        # Plot the Model Precipitation Biases
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
        # Denote that Model Data was Not Availabel
        bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                     horizontalalignment='center',verticalalignment='bottom',color='red',\
                     fontweight='semibold',fontsize=plt_param.TITLEfs);
        #######################################################################
        #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
        bias_ax.set_ylim(-5,5); bias_ax.set_yticks(np.arange(-5,5.1,2.5));
        bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
        #######################################################################
        #insert data statistics as title
        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
        BT.set_position(plt_param.Bloc_left);
        ###################################################################
        bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
        bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
        ###################################################################
        pr_fig.subplots_adjust(left = 0.125,right = 0.905,bottom = 0.155, top = 0.925)
        ###################################################################
        bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
                        fontsize = plt_param.TITLEfs)
        bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                        verticalalignment = 'bottom',horizontalalignment = 'right',\
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
    
