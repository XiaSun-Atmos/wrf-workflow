#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Measurement Data Exists

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
@modified: dmurray
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

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
import os,time
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
##########################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
from atmo_moisture import calc_dewpoint_from_MR
from atmo_moisture import calc_rh_from_MR
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
# Define SFC Measurement Heights 
SFCpress_hgt = 1;
SFCtemp_hgt = 2;
SFCrh_hgt = 2;
# per Jenni Kyrouac -- For the met systems: winds = 10m, pressure = 1m, temp/rh = 2m.

def psl_TPMXmodel(MDLin,OUTinfo):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth_meas = os.path.join(OUTinfo[0],'no_meas',DATE);
    if not os.path.exists(IMGpth_meas):
        os.makedirs(IMGpth_meas)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext',DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_none = os.path.join(OUTinfo[0],'no_data',DATE);
    if not os.path.exists(IMGpth_none):
        os.makedirs(IMGpth_none)
    ###########################################################################Z
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=img.imread(plt_param.NOAAlogo)
    FILpfx='pslmet_'
    
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
    psl_TRPplot(MEAStrp_plt,HRRRtrp_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""

    ###########################################################################
    """ Unpack Modeled Temperature, Pressure, and Mixing Ratio Data
    ###############################################################
    MDLin :: 
        HRRRtpmx_plt = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_tpmx)
        RAPtpmx_plt = (ESRLrap_ini,ESRLrap_loc,ESRLrap_xy,ESRLrap_tpmx)
        #-> each embedded tuple has length of initialization files
            ###################################################################
            MDLloc = (SGPsite,SGPlon,SGPlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLtpmx = depends on interpolation choice (see below)
                ESRLtpmx_o = (To,Tsfc_o,Po,Psfc_o,MXrat_o,MXrat_sfc_o,Prate_o,Pwtr_o,Pbl_o,Cbh_o,RHsfc_o)
                ESRLtpmx_int = (Tint,Tsfc_int,Pint,Psfc_int,MXrat_int,MXrat_sfc_int,Prate_int,Pwtr_int,Pbl_int,Cbh_int,RHsfc_int)
                ###############################################################
                ESRLtpmx = (ESRLtpmx_o,ESRLtpmx_int)
    ########################################################################"""
    # Model Initialization Times #        
    MDLini = MDLin[0]; #outputs model initialization time for each ini hour
    MDLxy = MDLin[2]; #outputs forecast time and measurement height
    MDLtpmx = MDLin[3] #outputs ini hour trp data (expand in loop)
    ###########################################################################
    MDLloc = MDLin[1][0]; #just store first initialization hour
    MDLloc_names = MDLloc[0];
    MDLlon_sites = MDLloc[1]; MDLlat_sites = MDLloc[2]; # desired grid point
    ######################################################################
    MDLlon_grids = MDLloc[3]; MDLlat_grids = MDLloc[4];# model grid point referenced
    MDLelev_grids = MDLloc[5]
    MDLdist2sites = MDLloc[6]
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
    #MEAStitle = MEASname+' ('+str(MEASalt)+' m)'
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    ###########################################################################

    ###########################################################################
    """ Measured TRP Information Not Available -- So Do Not Unpack """
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
    
                """ Extract Relevant Model Information """
                MDLtpmx_ini = MDLtpmx[INIref][INTPoi]
                ###################################################################
                # Extract Suface Estimates of Temperature, Relative Humidity, and 
                # Pressure
                MDLtemp = MDLtpmx_ini[1]
                MDLpres = MDLtpmx_ini[3]
                MDLmxrat = MDLtpmx_ini[5]
                MDLprecip = MDLtpmx_ini[11]
                MDLrh = MDLtpmx_ini[10]
                # surface--not height-values were extracted
                ###################################################################
                # Extract from teh Relevant Grid Point
                MDLtemp = MDLtemp[:,Gref];
                MDLpres = MDLpres[:,Gref];
                MDLmxrat = MDLmxrat[:,Gref];
                MDLprecip = MDLprecip[:,Gref];
                MDLrh = MDLrh[:,Gref];
                MDLtd = calc_dewpoint_from_MR(MDLmxrat/1000,MDLpres)
                MDLlcl = 125*(MDLtemp-MDLtd)
                # calculate RH if missing from MDL (e.g., RAP)
                if (np.isnan(MDLrh).all()):
                    MDLrh = calc_rh_from_MR(MDLmxrat/1000,MDLpres,MDLtemp)
                    MDLrh[np.where(MDLrh > 100)] = 100
                
                """ Define Relevant Plotting Specifics 
                --> Specifics will vary depending on whether an extended forecast 
                --> is being produced
                ###################################################################
                Extended is Necessary if Forecast Exceeds Model-Specific Standard 
                Length
                """
                if len(MDLtime) > FXlen: #indicating forecast hours exceed standard forecast length
                    """ Define Bounds for Extended Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLtemp_ext = MDLtemp; MDLpres_ext = MDLpres; MDLmxrat_ext = MDLmxrat;
                    MDLlcl_ext = MDLlcl; MDLrh_ext = MDLrh; MDLprecip_ext = MDLprecip;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    TEMPmin_ext,TEMPmax_ext,TEMPtick_ext = atmo_ylim((MDLtemp_ext),plt_param.TEMPrnd_bse,plt_param.TICKint)
                    PRESmin_ext,PRESmax_ext,PREStick_ext = atmo_ylim((MDLpres_ext),plt_param.PRESrnd_bse,plt_param.TICKint)
                    MXRATmin_ext,MXRATmax_ext,MXRATtick_ext = atmo_ylim((MDLmxrat_ext),plt_param.MXRATrnd_bse,plt_param.TICKint)
                    #RHmin_ext,RHmax_ext,RHtick_ext = atmo_ylim((MDLrh_ext),plt_param.RHrnd_bse,plt_param.TICKint)
                    RHmin_ext = 0; RHmax_ext = 110; RHtick_ext = np.arange(RHmin_ext,RHmax_ext+1,20)
                    LCLmin_ext,LCLmax_ext,LCLtick_ext = atmo_ylim((MDLlcl_ext),plt_param.LCLrnd_bse,plt_param.TICKint)
                    LCLmin_ext = 0;LCLmax_ext = max(1500,LCLmax_ext); LCLtick_ext = np.arange(LCLmin_ext,LCLmax_ext+1,plt_param.LCLrnd_bse);
                    PRmin_ext,PRmax_ext,PRtick_ext = atmo_bound_lim((MDLprecip_ext),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
                    APCPmin_ext,APCPmax_ext,APCPtick_ext = atmo_bound_lim((MDLprecip_ext),\
                                 plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin_ext = MDLtime_ext[0] - HRsec;
                    Xmax_ext = (MDLtime_ext[0]+EXTlen*HRsec) + HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick_ext = np.arange(MDLtime_ext[0],MDLtime_ext[0]+EXTlen*HRsec+1e-5,plt_param.Xdel_ext)
                    Xtick_ext_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick_ext]
                    Time_range_dates_ext = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick_ext]
                    Date_start_label_ext=Time_range_dates_ext[0].strftime('%Y-%m-%d')
                    Date_end_label_ext=Time_range_dates_ext[-1].strftime('%Y-%m-%d')
                    
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLtemp = MDLtemp[0:FXlen]; MDLpres = MDLpres[0:FXlen];
                    MDLmxrat = MDLmxrat[0:FXlen]; MDLrh = MDLrh[0:FXlen]; MDLlcl = MDLlcl[0:FXlen];
                    MDLprecip = MDLprecip[0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MDLtemp),plt_param.TEMPrnd_bse,plt_param.TICKint)
                    PRESmin,PRESmax,PREStick = atmo_ylim((MDLpres),plt_param.PRESrnd_bse,plt_param.TICKint)
                    MXRATmin,MXRATmax,MXRATtick = atmo_ylim((MDLmxrat),plt_param.MXRATrnd_bse,plt_param.TICKint)
                    #RHmin,RHmax,RHtick = atmo_ylim((MDLrh),plt_param.RHrnd_bse,plt_param.TICKint)
                    RHmin = 0; RHmax = 110; RHtick = np.arange(RHmin,RHmax+1,20)
                    LCLmin,LCLmax,LCLtick = atmo_ylim((MDLlcl),plt_param.LCLrnd_bse,plt_param.TICKint)
                    LCLmin = 0;LCLmax = max(1500,LCLmax); LCLtick = np.arange(LCLmin,LCLmax+1,plt_param.LCLrnd_bse);
                    PRmin,PRmax,PRtick = atmo_bound_lim((MDLprecip),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
                    APCPmin,APCPmax,APCPtick = atmo_bound_lim((MDLprecip),\
                                 plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0] - HRsec;
                    Xmax = (MDLtime[0]+(FXlen - 1)*HRsec) + HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                else:
                    """ Define Bounds for the Standard Forecast Period (i.e. + (FXlen - 1) Hrs) """
                    ###############################################################
                    # no model adjustment needed -- either equal to + (FXlen - 1) Fxcst or
                    # or less due to a download issue (i.e. MDLws = MDLws[0:FXlen])
                    ###############################################################
                    TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MDLtemp),plt_param.TEMPrnd_bse,plt_param.TICKint)
                    PRESmin,PRESmax,PREStick = atmo_ylim((MDLpres),plt_param.PRESrnd_bse,plt_param.TICKint)
                    MXRATmin,MXRATmax,MXRATtick = atmo_ylim((MDLmxrat),plt_param.MXRATrnd_bse,plt_param.TICKint)
                    #RHmin,RHmax,RHtick = atmo_ylim((MDLrh),plt_param.RHrnd_bse,plt_param.TICKint)
                    RHmin = 0; RHmax = 110; RHtick = np.arange(RHmin,RHmax+1,20)
                    LCLmin,LCLmax,LCLtick = atmo_ylim((MDLlcl),plt_param.LCLrnd_bse,plt_param.TICKint)
                    LCLmin = 0;LCLmax = max(1500,LCLmax); LCLtick = np.arange(LCLmin,LCLmax+1,plt_param.LCLrnd_bse);
                    ###############################################################
                    PRmin,PRmax,PRtick = atmo_bound_lim((MDLprecip),\
                                 plt_param.PRrnd_bse,plt_param.PRmin,plt_param.PRmax,\
                                 plt_param.PRrnd_bse)
                    APCPmin,APCPmax,APCPtick = atmo_bound_lim((MDLprecip),\
                                 plt_param.APCPrnd_bse,plt_param.APCPmin,plt_param.APCPmax,\
                                 plt_param.APCPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0] - HRsec;
                    Xmax = (MDLtime[0]+(FXlen - 1)*HRsec) + HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
    
                """  Model Data are Available -- Measurement Data Not Available """
                # would not enter loop if initialization time was not available
                """ ###############################################################
                Temperature Profile (instruct to share x-axis)
                ############################################################### """
                temp_fig, (temp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Temperature Time Histories
                #####################################
                temp_ax.plot(MDLtime,MDLtime*np.nan,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                temp_ax.plot(MDLtime,MDLtemp,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                temp_ax.plot(MDLtime,MDLtime*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ####################################################################
                handles, labels = temp_ax.get_legend_handles_labels()
                order = [1,0,2]
                temp_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                # Use Function to Denote Proper Temperature Bounds
                #TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MDLtemp),plt_param.TEMPrnd_bse,plt_param.TICKint)
                temp_ax.set_ylim(TEMPmin,TEMPmax); temp_ax.set_yticks(TEMPtick);
                temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                temp_ax.set_ylabel('Temperature, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});
                temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = temp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                """ ###############################################################
                # Plot the Model Temperature Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
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
                ####################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                temp_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                temp_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                            #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                imax = temp_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                temp_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
    
                """ ###############################################################
                #Pressure Plots (instruct to share x-axis)
                ############################################################### """
                pres_fig, (pres_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Atmospheric Pressure Time Histories
                ##############################################
                pres_ax.plot(MDLtime,MDLtime*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                pres_ax.plot(MDLtime,MDLtime*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                pres_ax.plot(MDLtime,MDLpres,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                pres_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###########################################################################
                handles, labels = pres_ax.get_legend_handles_labels()
                order = [1,0,2]
                pres_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###########################################################################
                # Use Function to Denote Proper Temperature Bounds
                #PRESmin,PRESmax,PREStick = atmo_ylim((MDLpres),plt_param.PRESrnd_bse,plt_param.TICKint)
                pres_ax.set_ylim(PRESmin,PRESmax); pres_ax.set_yticks(PREStick);
                pres_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###########################################################################
                pres_ax.set_ylabel('Pressure, mb',{'fontsize': plt_param.LABELfs});
                pres_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###########################################################################
                #insert data date as title
                T = pres_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Pressure Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
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
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                pres_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                pres_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                imax = pres_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_PRES_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pres_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
                
                """ ###############################################################
                #Lifting Condensation Level (LCL) Plots (instruct to share x-axis)
                ############################################################### """
                lcl_fig, (lcl_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Atmospheric LCL Time Histories
                ##############################################
                lcl_ax.plot(MDLtime,MDLtime*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                lcl_ax.plot(MDLtime,MDLtime*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                lcl_ax.plot(MDLtime,MDLlcl,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                lcl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###########################################################################
                handles, labels = lcl_ax.get_legend_handles_labels()
                order = [1,0,2]
                lcl_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###########################################################################
                # Use Function to Denote Proper LCL Bounds
                #LCLmin,LCLmax,LCLtick = atmo_ylim((MDLlcl),plt_param.LCLrnd_bse,plt_param.TICKint)
                #LCLmin = 0;LCLmax = max(1500,LCLmax_ext); LCLtick = np.arange(LCLmin,LCLmax+1,plt_param.LCLrnd_bse);
                lcl_ax.set_ylim(LCLmin,LCLmax); lcl_ax.set_yticks(LCLtick);
                lcl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###########################################################################
                lcl_ax.set_ylabel('Lifting Condensation \nLevel, m',{'fontsize': plt_param.LABELfs});
                lcl_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###########################################################################
                #insert data date as title
                T = lcl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model LCL Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
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
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-250,250); bias_ax.set_yticks(np.arange(-250,251,250));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, m',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                lcl_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                imax = lcl_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_LCL_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                lcl_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
                
                """ ###############################################################
                Mixing Ratio Profile (instruct to share x-axis)
                ############################################################### """
                mxrat_fig, (mxrat_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Mixing Ratio Time Histories
                #####################################
                mxrat_ax.plot(MDLtime,MDLtime*np.nan,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                mxrat_ax.plot(MDLtime,MDLmxrat,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                mxrat_ax.plot(MDLtime,MDLtime*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                mxrat_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ####################################################################
                handles, labels = mxrat_ax.get_legend_handles_labels()
                order = [1,0,2]
                mxrat_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                # Use Function to Denote Proper Mixing Ratio Bounds
                mxrat_ax.set_ylim(MXRATmin,MXRATmax); mxrat_ax.set_yticks(MXRATtick);
                mxrat_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                mxrat_ax.set_ylabel('Mixing Ratio, g kg-1',{'fontsize': plt_param.LABELfs});
                mxrat_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = mxrat_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                """ ###############################################################
                # Plot the Model Mixing Ratio Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
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
                ####################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                mxrat_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                mxrat_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                            #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                imax = mxrat_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_MXRAT_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                mxrat_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
    
                """ ###############################################################
                #Relative Humidity Plots (instruct to share x-axis)
                ############################################################### """
                rh_fig, (rh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Atmospheric Relative Humidity Time Histories
                ##############################################
                rh_ax.plot(MDLtime,MDLtime*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                rh_ax.plot(MDLtime,MDLtime*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                rh_ax.plot(MDLtime,MDLrh,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                rh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###########################################################################
                handles, labels = rh_ax.get_legend_handles_labels()
                order = [1,0,2]
                rh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###########################################################################
                # Use Function to Denote Proper Relative Humidity Bounds
                rh_ax.set_ylim(RHmin,RHmax); rh_ax.set_yticks(RHtick);
                rh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###########################################################################
                rh_ax.set_ylabel('Relative Humidity, %',{'fontsize': plt_param.LABELfs});
                rh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###########################################################################
                #insert data date as title
                T = rh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Relative Humidity Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)     
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
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-20,20); bias_ax.set_yticks(np.arange(-20,20.1,10));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                rh_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                rh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                imax = rh_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_RH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                rh_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
                
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
                pr_ax.bar(MDLtime,MDLtime*np.nan,width=120,color='blue',label='Obs',align='edge',linewidth=0.5)
                
                apcp_ax.plot(MDLtime,MDLtime*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='1-Hour Total')
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
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
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
                # Denote that Model Data was Not Availabel
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ####################################################################
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
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_PR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pr_fig.savefig(FILout,dpi=300); plt.close('all');

                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    """ Extended Temperature Forecast """
                    temp_fig, (temp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Temperature Time Histories
                    #####################################
                    #Observations 
                    temp_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                    temp_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    temp_ax.plot(MDLtime_ext,MDLtemp_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = temp_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    temp_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    temp_ax.set_ylim(TEMPmin_ext,TEMPmax_ext); temp_ax.set_yticks(TEMPtick_ext);
                    temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    temp_ax.set_ylabel('Temperature, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});
                    temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert data date as title
                    T = temp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                                 
                    """ ###########################################################
                    # Plot the Model Temperature Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #Define Bounds of Temperature - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    temp_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    temp_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                    imax = temp_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    temp_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
                    
                    """ Extended Pressure Forecasts """
                    pres_fig, (pres_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Atmospheric Pressure Time Histories
                    ##############################################
                    #Observations 
                    pres_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                    pres_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    pres_ax.plot(MDLtime_ext,MDLpres_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    pres_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = pres_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    pres_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    pres_ax.set_ylim(PRESmin_ext,PRESmax_ext); pres_ax.set_yticks(PREStick_ext);
                    pres_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    pres_ax.set_ylabel('Pressure, mb',{'fontsize': plt_param.LABELfs});
                    pres_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = pres_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Model Pressure Biases
                    ##########################################################  """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #Define Bounds of Pressure - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    pres_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    pres_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                    imax = pres_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_PRES_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    pres_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    """ Extended Lifting Condensation Level (LCL) Forecasts """
                    lcl_fig, (lcl_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the LCL Time Histories
                    ##############################################
                    #Observations 
                    lcl_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                    lcl_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    lcl_ax.plot(MDLtime_ext,MDLlcl_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    lcl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = lcl_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    lcl_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    lcl_ax.set_ylim(PRESmin_ext,PRESmax_ext); lcl_ax.set_yticks(PREStick_ext);
                    lcl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    lcl_ax.set_ylabel('Lifting Condensation \nLevel, m',{'fontsize': plt_param.LABELfs});
                    lcl_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = lcl_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Model Pressure Biases
                    ##########################################################  """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #Define Bounds of Pressure - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    lcl_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                    imax = lcl_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_LCL_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    lcl_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    """ Extended Mixing Ratio Forecast """
                    mxrat_fig, (mxrat_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Mixing Ratio Time Histories
                    #####################################
                    #Observations 
                    mxrat_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                    mxrat_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    mxrat_ax.plot(MDLtime_ext,MDLmxrat_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    mxrat_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = mxrat_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    mxrat_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    mxrat_ax.set_ylim(MXRATmin_ext,MXRATmax_ext); mxrat_ax.set_yticks(MXRATtick_ext);
                    mxrat_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    mxrat_ax.set_ylabel('Mixing Ratio, g kg-1',{'fontsize': plt_param.LABELfs});
                    mxrat_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert data date as title
                    T = mxrat_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                                 
                    """ ###########################################################
                    # Plot the Model Mixing Ratio Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #Define Bounds of Mixing Ratio - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, $\mathrm{deg}^{\circ}$ C',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    mxrat_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    mxrat_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                    imax = mxrat_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_MXRAT_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    mxrat_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
                    
                    """ Extended Relative Humidity Forecasts """
                    rh_fig, (rh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Atmospheric Relative Humidity Time Histories
                    ##############################################
                    #Observations 
                    rh_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                    rh_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    rh_ax.plot(MDLtime_ext,MDLrh_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    rh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = rh_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    rh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    rh_ax.set_ylim(RHmin_ext,RHmax_ext); rh_ax.set_yticks(RHtick_ext);
                    rh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    rh_ax.set_ylabel('Relative Humidity, %',{'fontsize': plt_param.LABELfs});
                    rh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = rh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Model Relative Humidity Biases
                    ##########################################################  """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #Define Bounds of Relative Humidity - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-20,20); bias_ax.set_yticks(np.arange(-20,20.1,10));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, mb',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    rh_ax.annotate('2 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    rh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                    imax = rh_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_RH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    rh_fig.savefig(FILout,dpi=300); plt.close('all');

                    """ Extended Precipitation Forecast """
                    pr_fig, (apcp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Precipitation Time Histories
                    #####################################
                    #Observations -- Only Plotting the QC'd Values
                    pr_ax=apcp_ax.twinx();
                    pr_ax.bar(MDLtime_ext,MDLtime_ext*np.nan,width=120,color='blue',label='Obs',align='edge',linewidth=0.5)
                    
                    apcp_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='1-Hour Total')
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
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #Define Bounds of Relative Humidity - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-5,5); bias_ax.set_yticks(np.arange(-5,5.1,2.5));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
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
