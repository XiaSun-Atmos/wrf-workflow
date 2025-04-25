#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case where both Measurement and Model Data are 
Available

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Precipitation
Plots Based on Measurements from the Disdrometers Installed at Various Sites
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
from datetime import datetime,timedelta
import os,time,warnings
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_beg import atmo_meas2hr_beg
from atmo_meas_to_hour_end import atmo_meas2hr_end
from atmo_meas_to_hour_upto import atmo_meas2hr_upto
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
##########################################
from atmo_bias import atmo_bias
from atmo_ylim_bias import atmo_blim
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
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
###############################################################################

def psl_PRplot(MEASin,MDLin,OUTinfo):
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
    FILpfx='psldisd_'
    
    """ Based on MDLname -- Define Relevant Variables """
    if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
        FXlen = fxcst_param.hrrr_v4_fx;
        EXTlen = fxcst_param.hrrr_v4_ext;
    elif MDLname == 'RAP_v5': # 'ESRL_RAP';
        FXlen = fxcst_param.rap_v5_fx;
        EXTlen = fxcst_param.rap_v5_ext;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    psl_TRPplot(MEAStrp_plt,HRRRtrp_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""

    ###########################################################################
    """ Unpack Modeled Temperature, Pressure, Mixing Ratio Data, and Precipitation
    ###############################################################
    MDLin :: 
        HRRRpr_plt = (PSLhrrr_ini,PSLhrrr_loc,PSLhrrr_xy,PSLhrrr_pr)
        RAPpr_plt = (PSLrap_ini,PSLrap_loc,PSLrap_xy,PSLrap_pr)
        #-> each embedded tuple has length of initialization files
            ###################################################################
            MDLloc = (PSLsite,PSLlon,PSLlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLpr = depends on interpolation choice (see below)
                PSLpr_o = (apcp_o)
                ###############################################################
                PSLpr = (PSLpr_o)
    ########################################################################"""
    # Model Initialization Times #        
    MDLini = MDLin[0]; #outputs model initialization time for each ini hour
    MDLxy = MDLin[2]; #outputs forecast time and measurement height
    MDLpr = MDLin[3] #outputs ini hour pr data (expand in loop)
    ######################################################################
    MDLloc = MDLin[1][0]; #just store first initialization hour
    ### Since we only have one station, we'll make arrays
    MDLloc_names = np.array([MDLloc[0]]);
    MDLlon_sites = np.array([MDLloc[1]]); MDLlat_sites = np.array([MDLloc[2]]); # desired grid point
    ######################################################################
    MDLlon_grids = np.array([MDLloc[3]]); MDLlat_grids = np.array([MDLloc[4]]);# model grid point referenced
    MDLelev_grids = np.array([MDLloc[5]]);
    MDLdist2sites = np.array([MDLloc[6]]);
    ###########################################################################
        
    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
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
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    ###########################################################################
    
    ###########################################################################
    """ Unpack Measured PR Information 
    ####################################
    MEASpr_plt = (MEASloc,MEAStme,MEASpr_in)       
    ---> MEASloc = (MEASlon,MEASlat,MEASalt)
    ---> MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus)
    ---> MEASpr_in = (MEASpr_minus,MEASpr,MEASpr_plus)
    ###########################################################################
    MEASpr = (MEASprate,MEASpr,MEASprsum)
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
    MEASprecip = atmo_cat(MEASpr_in[0],MEASpr_in[1],MEASpr_in[2],MEASpr_in[3],PRind_cat) 
    
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
            # Determine Whether Model Data Exists for Said Date
            if INIpos[ini] in MDLini:
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
                MDLpr_ini = MDLpr[INIref]
                ###################################################################
                # Extract Estimates of Precipitation
                MDLprecip = MDLpr_ini
                ###################################################################
                """ Define Relevant Plotting Specifics 
                --> Specifics will vary depending on whether an extended forecast 
                --> is being produced
                ###################################################################
                Extended is Necessary if Forecast Exceeds Model-Specific Standard 
                Length
                """
                if len(MDLtime) > FXlen: #indicating forecast hours exceed standard forecast length
                    """ Define Bounds for Extended Forecast Period """
                    ext_meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[0] + EXTlen*3600) 
                    ###############################################################
                    # Adjust Model Variables #
                    MDLprecip_ext = MDLprecip; 
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASprecip_ext_mu = atmo_meas2hr_upto(MEASprecip,None,MEAStime,MDLtime_ext,func='sum');
                    ###############################################################
                    PRmin_ext,PRmax_ext,PRtick_ext = atmo_bound_lim((MEASprecip[ext_meas_oi]),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
                    APCPmin_ext,APCPmax_ext,APCPtick_ext = atmo_bound_lim((MEASprecip_ext_mu,\
                                 MDLprecip_ext),plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin_ext = MDLtime_ext[0] - 3600;
                    Xmax_ext = (MDLtime_ext[0]+EXTlen*3600) + 3600;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick_ext = np.arange(MDLtime_ext[0],MDLtime_ext[0]+EXTlen*HRsec+1e-5,plt_param.Xdel_ext)
                    Xtick_ext_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick_ext]
                    Time_range_dates_ext = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick_ext]
                    Date_start_label_ext=Time_range_dates_ext[0].strftime('%Y-%m-%d')
                    Date_end_label_ext=Time_range_dates_ext[-1].strftime('%Y-%m-%d')
                    
                    """ Define Bounds for the Model-Specific Standard Forecast Period"""
                    meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[FXlen-1]) 
                    ###############################################################
                    # Adjust Model Variables #
                    MDLprecip = MDLprecip[0:FXlen]; 
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASprecip_mu = atmo_meas2hr_upto(MEASprecip,None,MEAStime,MDLtime,func='sum');
                    ###############################################################
                    PRmin,PRmax,PRtick = atmo_bound_lim((MEASprecip[meas_oi]),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
                    APCPmin,APCPmax,APCPtick = atmo_bound_lim((MEASprecip_mu,\
                                 MDLprecip),plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0] - 3600;
                    Xmax = (MDLtime[0]+(FXlen - 1)*3600) + 3600;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                else:
                    """ Define Bounds for the Model-Specific Standard Forecast Period"""
                    ###############################################################
                    # Know that an Extended Forecast is Non-Existent. Therefore,
                    # Simply Define meas_oi based on the Desired Forecast Time of 
                    # 86400
                    ###############################################################
                    meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[0]+(FXlen - 1)*3600)
                    ###############################################################
                    # no model adjustment needed -- either equal to + (FXlen - 1) Fxcst or
                    # or less due to a download issue (i.e. MDLws = MDLws[0:FXlen])
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASprecip_mu = atmo_meas2hr_upto(MEASprecip,None,MEAStime,MDLtime,func='sum');
                    ###############################################################
                    PRmin,PRmax,PRtick = atmo_bound_lim((MEASprecip[meas_oi]),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
                    APCPmin,APCPmax,APCPtick = atmo_bound_lim((MEASprecip_mu),\
                                 plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0] - 3600;
                    Xmax = (MDLtime[0]+(FXlen - 1)*3600) + 3600;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
    
                ###################################################################
                """ Measurements and Model Data are Available """
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
                apcp_ax.plot(MDLtime,MDLprecip,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
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
                APCPbias = MDLprecip - MEASprecip_mu
                bias_ax.plot(MDLtime,APCPbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                #Define Bounds of Precipitation - Model Bias to Determine Y-Bounds
                APCPb_min,APCPb_max,APCPb_tick = atmo_bound_lim((APCPbias),\
                                 plt_param.APCPbias_bse,plt_param.APCPbmin,plt_param.APCPbmax,\
                                 plt_param.APCPbias_bse)
                bias_ax.set_ylim(APCPb_min,APCPb_max); bias_ax.set_yticks(APCPb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                RMSEpr,BIASpr = atmo_bias(APCPbias)
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: ' + str(RMSEpr)+' mm | Mean Bias: '+\
                                      str(BIASpr)+' mm',{'fontsize': plt_param.TITLEfs,\
                                      'horizontalalignment':'left'},pad=2.5) 
                BT.set_position(plt_param.Bloc_left);
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
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = pr_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_PR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pr_fig.savefig(FILout,dpi=300); plt.close('all');
    
                
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    """ Extended Precipitation Forecast """
                    pr_fig, (apcp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Precipitation Time Histories
                    #####################################
                    #Observations -- Only Plotting the QC'd Values
                    pr_ax=apcp_ax.twinx();
                    pr_ax.bar(MEAStime,MEASprecip,width=120,color='blue',label='Obs',align='edge',linewidth=0.5)
                    
                    apcp_ax.plot(MDLtime_ext,MEASprecip_ext_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='1-Hour Total')
                    apcp_ax.plot(MDLtime_ext,MDLprecip_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    apcp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    prlines,prlabels = pr_ax.get_legend_handles_labels()
                    apcplines,apcplabels = apcp_ax.get_legend_handles_labels()
                    apcp_ax.legend(prlines+apcplines, prlabels+apcplabels, loc=1,\
                                   fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    apcp_ax.set_ylim(APCPmin_ext,APCPmax_ext); apcp_ax.set_yticks(APCPtick_ext);
                    apcp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    pr_ax.set_ylim(PRmin_ext,PRmax_ext); pr_ax.set_yticks(PRtick_ext);
                    pr_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    apcp_ax.set_ylabel('1-Hr Total Precip, mm',{'fontsize': plt_param.LABELfs});
                    apcp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    pr_ax.set_ylabel('2-min Obs Precip, mm',{'fontsize': plt_param.LABELfs});
                    ###############################################################
                    #insert station info as title
                    T = apcp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                                 
                    """ ###########################################################
                    # Plot the Model Precip Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    APCPbias = MDLprecip_ext - MEASprecip_ext_mu
                    bias_ax.plot(MDLtime_ext,APCPbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                 ncol = 1,frameon = False);    
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                    ###############################################################
                    #Define Bounds of Temperature - Model Bias to Determine Y-Bounds
                    #APCPb_min,APCPb_max,APCPb_tick = atmo_blim((APCPbias),plt_param.APCPbias_bse,plt_param.TICKint)
                    APCPb_min,APCPb_max,APCPb_tick = atmo_bound_lim((APCPbias),\
                                 plt_param.APCPbias_bse,plt_param.APCPbmin,plt_param.APCPbmax,\
                                 plt_param.APCPbias_bse)
                    bias_ax.set_ylim(APCPb_min,APCPb_max); bias_ax.set_yticks(APCPb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    RMSEpr,BIASpr = atmo_bias(APCPbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEpr)+' mm | Mean Bias: '+\
                                          str(BIASpr)+' mm',{'fontsize': plt_param.TITLEfs,\
                                          'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    pr_fig.subplots_adjust(left = 0.125,right = 0.905,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = pr_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_PR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    pr_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
                    # Delete Extended Information to Preclude Incorrect Plotting
                    del MDLtime_ext
                    
            else:
                """ Initialization Time Not Available """
                ###################################################################
                MDLpseudo_time = np.arange(0,(FXlen - 1)*3600+1e-5,step=3600)
                MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
                ###################################################################
                #Define Xtick Information (Only Needs to be Defined Once)
                Xmin = MDLtime[0] - 3600;
                Xmax = (MDLtime[0]+(FXlen - 1)*3600) + 3600;
                ###################################################################
                #Define Xtick Information (Only Needs to be Defined Once)
                Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
                Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                
                ###################################################################
                """ Perform Model Averaging for Plotting and Statistical Analyses """
                MEASprecip_mu = atmo_meas2hr_upto(MEASprecip,None,MEAStime,MDLtime,func='sum');
                ###################################################################
                """ Define Relevant Plotting Specifics """
                # Define Proper Indices to Determine Axes Parameters
                meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
                ###################################################################
                PRmin,PRmax,PRtick = atmo_bound_lim((MEASprecip[meas_oi]),\
                             plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                             plt_param.PRrnd_bse)
                APCPmin,APCPmax,APCPtick = atmo_bound_lim((MEASprecip_mu),\
                             plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                             plt_param.APCPrnd_bse)
                """ ###############################################################
                Precipitation Profile (instruct to share x-axis)
                ############################################################### """
                pr_fig, (apcp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ####################################################################
                # Plot the Precipitation Time Histories
                #####################################
                #Observations -- Only Plotting the QC'd Values
                pr_ax = apcp_ax.twinx()
                pr_ax.bar(MEAStime,MEASprecip,width=120,color='blue',label='Obs',align='edge',linewidth=0.5)
                
                apcp_ax.plot(MDLtime,MEASprecip_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Total')
                apcp_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
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
                # Plot the Model Precip Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
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
                ###################################################################
                # Denote that Model Data was Not Availabel
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-5,5); bias_ax.set_yticks(np.arange(-5,5.1,2.5));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
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
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
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
                ###################################################################
    
