#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Measurement Data Exists
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
  neighbouring NaNs appear with some weird black structure indicating that the 
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
import matplotlib.gridspec as gridspec
#####################################
from coord_distance import dist2coord
#####################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
from atmo_spatial_ylim import atmo_sptl_ylim
from wd_spatial_ylim import wd_sptl_ylim
###################################################
from wd_to_cmap import wd_cmap
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

###############################################################################
# Define Desired Height Limit
MINlow = 0; MAXlow = 2950;
#MINhgh = 1250; MAXhgh = 4250;
MINhgh = 0; MAXhgh = 6097;

def rwp_Wmodel(MDLin,OUTinfo,frequency=915):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth_meas = os.path.join(OUTinfo[0], 'no_meas',DATE);
    if not os.path.exists(IMGpth_meas):
        os.makedirs(IMGpth_meas)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext', DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_none = os.path.join(OUTinfo[0],'no_data', DATE);
    if not os.path.exists(IMGpth_none):
        os.makedirs(IMGpth_none)
    ##################################################
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=image.imread(plt_param.NOAAlogo);
    FILpfx='rwpwinds'+str(frequency)+'_'
    
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
    """ Measured Wind Information Not Available -- Do Not Unpack 
    --> Measurement Heights Required for Proper Barb Information (Defined Below)
    ####################################################################### """
    # Define Measurement Arrays for Proper Reference in Plots
    MEAShgt_low = np.asarray([125.4638,146.6121,167.7604,188.9087,210.057,231.2053,\
                    252.3535,273.5018,294.6501,315.7984,336.9467,358.095,379.2433,\
                    400.3915,421.5398,442.6881,463.8364,484.9847,506.133,527.2812,\
                    548.4295,569.5778,590.726,611.8743,633.0225,654.1708,675.319,\
                    696.4673,717.6155,738.7638,759.912,781.0603,802.2086,823.3568,\
                    844.5051,865.6533,886.8016,907.9498,929.0981,950.2463,971.3946,\
                    992.5428,1013.6911,1034.8394,1055.9877,1077.136,1098.2843,1119.4326,\
                    1140.5809,1161.7292,1182.8776,1204.0259,1225.1742,1246.3225,1267.4708,\
                    1288.6191,1309.7675,1330.9158,1352.0641,1373.2124,1394.3607,1415.509,\
                    1436.6573,1457.8057,1478.954,1500.1023,1521.2506,1542.3989,1563.5472,\
                    1584.6956,1605.8439,1626.9922,1648.1405,1669.2888,1690.4371,1711.5854,\
                    1732.7338,1753.8821,1775.0304,1796.1787,1817.327,1838.4753,1859.6237,1880.772,\
                    1901.9203,1923.0686,1944.2169,1965.3652,1986.5135,2007.6619,2028.8102,2049.9585,\
                    2071.1067,2092.255,2113.403,2134.5513,2155.6995,2176.8477,2197.9958,2219.144,\
                    2240.2922,2261.4404,2282.5886,2303.7368,2324.885,2346.0332,2367.1814,2388.3296,\
                    2409.4778,2430.626,2451.7742,2472.9224,2494.0706,2515.2188,2536.367,2557.5151,\
                    2578.6633,2599.8115,2620.9597,2642.108,2663.256])
    MEAShgt_hgh = np.asarray([1511.2052,1582.052,1652.8988,1723.7456,1794.5924,1865.4392,1936.286,\
                    2007.1328,2077.9795,2148.8262,2219.6729,2290.5195,2361.3662,2432.213,2503.0596,\
                    2573.9062,2644.753,2715.5996,2786.4463,2857.293,2928.1396,2998.9863,3069.833,\
                    3140.6797,3211.5264,3282.373,3353.2197,3424.0664,3494.913,3565.7598,3636.6064,\
                    3707.4531,3778.2998,3849.1465,3919.9932,3990.8398,4061.6865,4132.533,4203.38,\
                    4274.2266,4345.073,4415.92,4486.7666,4557.6133,4628.46])
    
    ###########################################################################
    """ Define Reference Indices for Barb Visualization """
    MEAShgt_low_max = MEAShgt_low.shape[0];
    MEAShgt_hgh_max = MEAShgt_hgh.shape[0];
    MEAShgt_low_min = np.where(MEAShgt_low >= MINlow)[0][0] 
    MEAShgt_hgh_min = np.where(MEAShgt_hgh >= MINhgh)[0][0]
    
    #######################################################################
    # set up the colormaps
    cmap_ws = get_var_cmap('ws')
    cmap_bias = get_var_cmap('bias')
    cmap_wd = get_var_cmap('wd')
    cmap_snr = get_var_cmap('snr')

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
                ############################################################### """
                # Maximum and Minimum Measurement Heights 
                HGTmin_low,HGTmax_low,HGTtick_low = atmo_ylim((np.asarray([MINlow,MAXlow])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTmin_hgh,HGTmax_hgh,HGTtick_hgh = atmo_ylim((np.asarray([MINhgh,MAXhgh])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTtick_low_label = [str(ht/1000) for ht in HGTtick_low]
                HGTtick_hgh_label = [str(ht/1000) for ht in HGTtick_hgh]
                #######################################################################
                # Define the Max Height Ind
                #MEAShgt_low_max = np.where(MEAShgt_low > MAXlow)[0][0] 
                #MEAShgt_hgh_max = np.where(MEAShgt_hgh > MAXhgh)[0][0]
                ###################################################################
                # Measurement Barbs Need Not be Defined
                
                if len(MDLtime) > FXlen:
                    """ Define Bounds for Extended Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLws_ext = MDLws; MDLwd_ext = MDLwd;
                    MDLu_ext = MDLu; MDLv_ext = MDLv;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    # Define Relevant Model Grid
                    MDLtme_ext_mesh, MDLhgt_ext_mesh = np.meshgrid(MDLtime_ext,MDLhgt); #define reference meshgrid
                    MDLgrid_ext = (MDLtime_ext[0],MDLtime_ext[-1],MDLhgt[0],MDLhgt[-1])
                    ###############################################################
                    low_ext_mdl_oi = np.logical_and(MDLhgt_ext_mesh >= MINlow, MDLhgt_ext_mesh <= MAXlow)
                    hgh_ext_mdl_oi = np.logical_and(MDLhgt_ext_mesh >= MINhgh, MDLhgt_ext_mesh <= MAXhgh)
                    ###############################################################
                    # pass through array indices only of interest for y-lims
                    #WSmin_low_ext,WSmax_low_ext,WStick_low_ext = atmo_sptl_ylim((MDLws_ext[low_ext_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    #WSmin_hgh_ext,WSmax_hgh_ext,WStick_hgh_ext = atmo_sptl_ylim((MDLws_ext[hgh_ext_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    WSmin_low_ext,WSmax_low_ext,WStick_low_ext = atmo_bound_lim((MDLws_ext[low_ext_mdl_oi]),
                                    plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,plt_param.WSrnd_bse)
                    WSmin_hgh_ext,WSmax_hgh_ext,WStick_hgh_ext = atmo_bound_lim((MDLws_ext[hgh_ext_mdl_oi]),
                                    plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,plt_param.WSrnd_bse)
                    
                    #WDmin_low_ext,WDmax_low_ext,WDtick_low_ext = wd_sptl_ylim((MDLwd_ext[low_ext_mdl_oi]));
                    #WDmin_hgh_ext,WDmax_hgh_ext,WDtick_hgh_ext = wd_sptl_ylim((MDLwd_ext[hgh_ext_mdl_oi]));   
                    WDmin_low_ext=0;WDmax_low_ext=360;WDtick_low_ext=np.arange(0,361,60);
                    WDmin_hgh_ext=0;WDmax_hgh_ext=360;WDtick_hgh_ext=np.arange(0,361,60);
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin_ext = MDLtime_ext[0];
                    Xmax_ext = (MDLtime_ext[0]+EXTlen*HRsec);
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
                    MDLws = MDLws[:,0:FXlen]; MDLwd = MDLwd[:,0:FXlen];
                    MDLu = MDLu[:,0:FXlen]; MDLv = MDLv[:,0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    # Define Relevant Model Grid
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    MDLgrid = (MDLtime[0],MDLtime[-1],MDLhgt[0],MDLhgt[-1])
                    ###############################################################
                    low_mdl_oi = np.logical_and(MDLhgt_mesh >= MINlow, MDLhgt_mesh <= MAXlow)
                    hgh_mdl_oi = np.logical_and(MDLhgt_mesh >= MINhgh, MDLhgt_mesh <= MAXhgh)
                    ###############################################################
                    # pass through array indices only of interest for y-lims
                    #WSmin_low,WSmax_low,WStick_low = atmo_sptl_ylim((MDLws[low_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    #WSmin_hgh,WSmax_hgh,WStick_hgh = atmo_sptl_ylim((MDLws[hgh_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    WSmin_low,WSmax_low,WStick_low = atmo_bound_lim((MDLws[low_mdl_oi]),plt_param.WSrnd_bse,
                                plt_param.WSmin,plt_param.WSmax,plt_param.WSrnd_bse)
                    WSmin_hgh,WSmax_hgh,WStick_hgh = atmo_bound_lim((MDLws[hgh_mdl_oi]),
                                plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,plt_param.WSrnd_bse)
                    
                    #WDmin_low,WDmax_low,WDtick_low = wd_sptl_ylim((MDLwd[low_mdl_oi]));
                    #WDmin_hgh,WDmax_hgh,WDtick_hgh = wd_sptl_ylim((MDLwd[hgh_mdl_oi]));   
                    WDmin_low=0;WDmax_low=360;WDtick_low=np.arange(0,361,60);
                    WDmin_hgh=0;WDmax_hgh=360;WDtick_hgh=np.arange(0,361,60);
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0];
                    Xmax = MDLtime[0]+(FXlen-1)*HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen-1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick_ext]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                    """ Define Model Barb Information (Low Measurement Height Mode)
                    ########################################################### """
                    BARBmdl_Loi = np.zeros(MDLtme_mesh.shape,dtype = bool)
                    BARBmdl_ext_Loi = np.zeros(MDLtme_ext_mesh.shape,dtype = bool)
                    ###############################################################
                    # define the measurement heights of interest
                    BARBhgt_Loi = MEAShgt_low[np.arange(MEAShgt_low_min,MEAShgt_low_max,plt_param.BARBrwp_Ly)]
                    ###############################################################
                    for hoi in BARBhgt_Loi:
                         # identify the measurement height reference indice #
                        barb_hgt = np.where(np.abs(MDLhgt - hoi) - np.min(np.abs(MDLhgt - hoi)) \
                                             < 1e-5)[0]
                        ###########################################################
                        # y-steps are the same whether measurement or model
                        for i in np.arange(plt_param.BARBrwp_x_base,len(MDLtime),plt_param.BARBrwp_x):
                            BARBmdl_Loi[barb_hgt,i] = True
                        ###############################################################
                        for i in np.arange(plt_param.BARBrwp_x_base,len(MDLtime_ext),plt_param.BARBrwp_x_ext):
                            BARBmdl_ext_Loi[barb_hgt,i] = True
                            
                    """ Define Model Barb Information (High Measurement Height Mode)
                    ########################################################### """
                    BARBmdl_Hoi = np.zeros(MDLtme_mesh.shape,dtype = bool)
                    BARBmdl_ext_Hoi = np.zeros(MDLtme_ext_mesh.shape,dtype = bool)
                    ###############################################################
                    # define the measurement heights of interest
                    BARBhgt_Hoi = MEAShgt_hgh[np.arange(MEAShgt_hgh_min,MEAShgt_hgh_max,plt_param.BARBrwp_Hy)]
                    ###############################################################
                    for hoi in BARBhgt_Hoi:
                         # identify the measurement height reference indice #
                        barb_hgt = np.where(np.abs(MDLhgt - hoi) - np.min(np.abs(MDLhgt - hoi)) \
                                             < 1e-5)[0]
                        ###########################################################
                        # y-steps are the same whether measurement or model
                        for i in np.arange(plt_param.BARBrwp_x_base,len(MDLtime),plt_param.BARBrwp_x):
                            BARBmdl_Hoi[barb_hgt,i] = True
                        ###############################################################
                        for i in np.arange(plt_param.BARBrwp_x_base,len(MDLtime_ext),plt_param.BARBrwp_x_ext):
                            BARBmdl_ext_Hoi[barb_hgt,i] = True
                    
                    """ Define Wind Direction Colormap Based on Array Range
                    ########################################################### 
                    cmap_wd_low = wd_cmap(WDmin_low,WDmax_low);
                    cmap_wd_hgh = wd_cmap(WDmin_hgh,WDmax_hgh);
                    ###############################################################
                    cmap_wd_low_ext = wd_cmap(WDmin_low_ext,WDmax_low_ext);
                    cmap_wd_hgh_ext = wd_cmap(WDmin_hgh_ext,WDmax_hgh_ext);
                    """
                    ########################################################### 
                    cmap_wd_low = cmap_wd
                    cmap_wd_hgh = cmap_wd
                    ###############################################################
                    cmap_wd_low_ext = cmap_wd
                    cmap_wd_hgh_ext = cmap_wd
                    
                else:
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # no model adjustment needed
                    # or less due to a download issue (i.e. MDLws = MDLws[0:FXlen])
                    ###############################################################
                    # Define Relevant Model Grid
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    MDLgrid = (MDLtime[0],MDLtime[-1],MDLhgt[0],MDLhgt[-1])
                    ###############################################################
                    low_mdl_oi = np.logical_and(MDLhgt_mesh >= MINlow, MDLhgt_mesh <= MAXlow)
                    hgh_mdl_oi = np.logical_and(MDLhgt_mesh >= MINhgh, MDLhgt_mesh <= MAXhgh)
                    ###############################################################
                    # pass through array indices only of interest for y-lims
                    #WSmin_low,WSmax_low,WStick_low = atmo_sptl_ylim((MDLws[low_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    #WSmin_hgh,WSmax_hgh,WStick_hgh = atmo_sptl_ylim((MDLws[hgh_mdl_oi]),plt_param.WSrnd_bse,plt_param.TICKint)    
                    WSmin_low,WSmax_low,WStick_low = atmo_bound_lim((MDLws[low_mdl_oi]),plt_param.WSrnd_bse,
                                plt_param.WSmin,plt_param.WSmax,plt_param.WSrnd_bse)
                    WSmin_hgh,WSmax_hgh,WStick_hgh = atmo_bound_lim((MDLws[hgh_mdl_oi]),
                                plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,plt_param.WSrnd_bse)
                    
                    #WDmin_low,WDmax_low,WDtick_low = wd_sptl_ylim((MDLwd[low_mdl_oi]));
                    #WDmin_hgh,WDmax_hgh,WDtick_hgh = wd_sptl_ylim((MDLwd[hgh_mdl_oi]));   
                    WDmin_low=0;WDmax_low=360;WDtick_low=np.arange(0,361,60);
                    WDmin_hgh=0;WDmax_hgh=360;WDtick_hgh=np.arange(0,361,60);
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0];
                    Xmax = MDLtime[0]+(FXlen-1)*HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen-1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')                
    
                    """ Define Model Barb Information
                    ########################################################### """
                    BARBmdl_Loi = np.zeros(MDLtme_mesh.shape,dtype = bool)
                    BARBmdl_Hoi = np.zeros(MDLtme_mesh.shape,dtype = bool)
                    ###############################################################
                    # define the measurement heights of interest
                    BARBhgt_Loi = MEAShgt_low[np.arange(MEAShgt_low_min,MEAShgt_low_max,plt_param.BARBrwp_Ly)]
                    BARBhgt_Hoi = MEAShgt_hgh[np.arange(MEAShgt_hgh_min,MEAShgt_hgh_max,plt_param.BARBrwp_Hy)]
                    ###############################################################
                    for i in np.arange(plt_param.BARBrwp_x_base,len(MDLtime),plt_param.BARBrwp_x):
                        for hoi in BARBhgt_Loi:
                             # identify the measurement height reference indice #
                            barb_hgt = np.where(np.abs(MDLhgt - hoi) - \
                                        np.min(np.abs(MDLhgt - hoi)) < 1e-5)[0]
                            #######################################################
                            BARBmdl_Loi[barb_hgt,i] = True
                        ###############################################################
                        for hoi in BARBhgt_Hoi:
                             # identify the measurement height reference indice #
                            barb_hgt = np.where(np.abs(MDLhgt - hoi) - \
                                                np.min(np.abs(MDLhgt - hoi)) < 1e-5)[0]
                            ###########################################################
                            BARBmdl_Hoi[barb_hgt,i] = True
                    
                    """ Define Wind Direction Colormap Based on Array Range
                    ########################################################### 
                    cmap_wd_low = wd_cmap(WDmin_low,WDmax_low);
                    cmap_wd_hgh = wd_cmap(WDmin_hgh,WDmax_hgh);
                    """
                    cmap_wd_low = cmap_wd
                    cmap_wd_hgh = cmap_wd
    
                """ ###############################################################
                Wind Speed Plots (High [i.e. Low-Resolution] Mode)
                ############################################################### """
                ws_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                # ws_fig, (rwp_ax,mdl_ax,bias_ax) = gridspec.GridSpec(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                #          height_ratios=[2,2,1],\
                #          facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Speed Time Histories
                ##########################################
                ws_fig.suptitle(str(frequency)+' MHz RWP WS',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    ws_rwp = rwp_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=WSmin_low,vmax=WSmax_low,\
                                         origin='lower',aspect='auto',cmap=cmap_ws);    
                rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = ws_fig.colorbar(ws_rwp,ax=rwp_ax,ticks = WStick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
                rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                #######################################################################
                # set figure y-label
                rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Measurement Data was Not Available
                rwp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #insert stn name/alt as title
                T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                ###################################################################
                # Plot the Model Wind Speed Time Histories
                ##########################################
                ws_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLws,vmin = WSmin_low,\
                                            vmax = WSmax_low,shading='gouraud',cmap=cmap_ws)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                mdl_ax.barbs(MDLtme_mesh[BARBmdl_Loi],MDLhgt_mesh[BARBmdl_Loi],MDLu[BARBmdl_Loi],\
                            MDLv[BARBmdl_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                mdl_ax.set_ylim(HGTmin_low,HGTmax_low); mdl_ax.set_yticks(HGTtick_low);mdl_ax.set_yticklabels(HGTtick_low_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WS'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ws = bias_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin_low,HGTmax_low); bias_ax.set_yticks(HGTtick_low);bias_ax.set_yticklabels(HGTtick_low_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                #######################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                #######################################################################
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
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
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_WSlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Wind Direction Plots (High [i.e. Low-Resolution] Mode)
                ############################################################### """
                wd_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Direction Time Histories
                ##############################################
                wd_fig.suptitle(str(frequency)+' MHz RWP WD',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                                fontsize=plt_param.TITLEfs);
                ###################################################################
                wd_rwp = rwp_ax.imshow(MDLwd*np.nan,extent=MDLgrid,vmin=WDmin_low,vmax=WDmax_low,\
                                     origin='lower',aspect='auto',cmap=cmap_wd_low);  
                rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = wd_fig.colorbar(wd_rwp,ax=rwp_ax,ticks = WDtick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
                rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                #######################################################################
                # set figure y-label
                rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                rwp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #insert stn name/alt as title
                T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Wind Direction Time Histories
                ##############################################
                wd_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLwd,vmin = WDmin_low,\
                                            vmax = WDmax_low,shading='gouraud',cmap=cmap_wd_low);
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_low,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                mdl_ax.barbs(MDLtme_mesh[BARBmdl_Loi],MDLhgt_mesh[BARBmdl_Loi],MDLu[BARBmdl_Loi],\
                            MDLv[BARBmdl_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                mdl_ax.set_ylim(HGTmin_low,HGTmax_low); mdl_ax.set_yticks(HGTtick_low);mdl_ax.set_yticklabels(HGTtick_low_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WD'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Direction Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_wd = bias_ax.imshow(MDLwd*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin_low,HGTmax_low); bias_ax.set_yticks(HGTtick_low);bias_ax.set_yticklabels(HGTtick_low_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                #######################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                #######################################################################
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
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
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_WDlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Signal-to-Noise Ratio Plots (High [i.e. Low-Resolution] Mode)
                ############################################################### """
                #snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[3,1,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Signal-to-Noise Time Histories
                ###############################################
                snr_fig.suptitle(str(frequency)+' MHz RWP Range-Corrected SNR',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    snr_rwp = rwp_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=-40,vmax=10,\
                                         origin='lower',aspect='auto',cmap=cmap_snr);   
                rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = snr_fig.colorbar(snr_rwp,ax=rwp_ax,ticks = np.arange(-40,10.1,10),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect3)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('SNR, dB', fontsize=plt_param.LABELfs);
                ###################################################################
                rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
                rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                #######################################################################
                # set figure y-label
                rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                rwp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #insert stn name/alt as title
                T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################          
                                     
                """ ###############################################################
                Plot the Model Wind Speed Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ws = ws_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_ws);   
                ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(bias_ws,ax=ws_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                ws_ax.set_ylim(HGTmin_low,HGTmax_low); ws_ax.set_yticks(HGTtick_low);ws_ax.set_yticklabels(HGTtick_low_label);
                ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                ws_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                ws_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = ws_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################                     
                                     
                """ ###############################################################
                Define and Plot the Model Wind Direction Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_wd = wd_ax.imshow(MDLwd*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_wd);   
                wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(bias_wd,ax=wd_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                ###################################################################
                wd_ax.set_ylim(HGTmin_low,HGTmax_low); wd_ax.set_yticks(HGTtick_low);wd_ax.set_yticklabels(HGTtick_low_label);
                wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                #######################################################################
                # set figure y-label
                wd_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                wd_ax.set_xticks(Xtick); wd_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    wd_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                wd_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                wd_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                wd_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = wd_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                snr_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                wd_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                wd_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = snr_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_SNRlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                snr_fig.savefig(FILout,dpi=300); plt.close('all');         
                ###################################################################
                ###################################################################
                ###################################################################
                """ ###############################################################
                Wind Speed Plots (Low [i.e. High-Resolution] Mode)
                ############################################################### """
                ws_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Speed Time Histories
                ##########################################
                ws_fig.suptitle(str(frequency)+' MHz RWP WS',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    ws_rwp = rwp_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=WSmin_hgh,vmax=WSmax_hgh,\
                                         origin='lower',aspect='auto',cmap=cmap_ws);    
                rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = ws_fig.colorbar(ws_rwp,ax=rwp_ax,ticks = WStick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
                rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                #######################################################################
                # set figure y-label
                rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Measurement Data was Not Available
                rwp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #insert stn name/alt as title
                T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                ###################################################################
                # Plot the Model Wind Speed Time Histories
                ##########################################
                ws_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLws,vmin = WSmin_hgh,\
                                            vmax = WSmax_hgh,shading='gouraud',cmap=cmap_ws)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                mdl_ax.barbs(MDLtme_mesh[BARBmdl_Hoi],MDLhgt_mesh[BARBmdl_Hoi],MDLu[BARBmdl_Hoi],\
                            MDLv[BARBmdl_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                mdl_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); mdl_ax.set_yticks(HGTtick_hgh);mdl_ax.set_yticklabels(HGTtick_hgh_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WS'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ws = bias_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); bias_ax.set_yticks(HGTtick_hgh);bias_ax.set_yticklabels(HGTtick_hgh_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                #######################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_WShgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Wind Direction Plots (Low [i.e. High-Resolution] Mode)
                ############################################################### """
                wd_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Wind Direction Time Histories
                ##############################################
                wd_fig.suptitle(str(frequency)+' MHz RWP WD',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                                fontsize=plt_param.TITLEfs);
                ###################################################################
                wd_rwp = rwp_ax.imshow(MDLwd*np.nan,extent=MDLgrid,vmin=WDmin_hgh,vmax=WDmax_hgh,\
                                     origin='lower',aspect='auto',cmap=cmap_wd_hgh);  
                rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = wd_fig.colorbar(wd_rwp,ax=rwp_ax,ticks = WDtick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
                rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                #######################################################################
                # set figure y-label
                rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                rwp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #insert stn name/alt as title
                T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Wind Direction Time Histories
                ##############################################
                wd_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLwd,vmin = WDmin_hgh,\
                                            vmax = WDmax_hgh,shading='gouraud',cmap=cmap_wd_hgh);    
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_hgh,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                ###################################################################
                # overlay barbs
                mdl_ax.barbs(MDLtme_mesh[BARBmdl_Hoi],MDLhgt_mesh[BARBmdl_Hoi],MDLu[BARBmdl_Hoi],\
                            MDLv[BARBmdl_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
                ###################################################################
                mdl_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); mdl_ax.set_yticks(HGTtick_hgh);mdl_ax.set_yticklabels(HGTtick_hgh_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' WD'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Wind Direction Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_wd = bias_ax.imshow(MDLwd*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); bias_ax.set_yticks(HGTtick_hgh);bias_ax.set_yticklabels(HGTtick_hgh_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                #######################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
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
                #######################################################################
                # add a logo
                imax = wd_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_WDhgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Signal-to-Noise Ratio Plots (Low [i.e. High-Resolution] Mode)
                ############################################################### """
                #snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[3,1,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Signal-to-Noise Time Histories
                ###############################################
                snr_fig.suptitle(str(frequency)+' MHz RWP Range-Corrected SNR',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    snr_rwp = rwp_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=-40,vmax=10,\
                                         origin='lower',aspect='auto',cmap=cmap_snr);   
                rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = snr_fig.colorbar(snr_rwp,ax=rwp_ax,ticks = np.arange(-40,10.1,10),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect3)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('SNR, dB', fontsize=plt_param.LABELfs);
                ###################################################################
                rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
                rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                #######################################################################
                # set figure y-label
                rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                rwp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #insert stn name/alt as title
                T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################          
                                     
                """ ###############################################################
                Plot the Model Wind Speed Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ws = ws_ax.imshow(MDLws*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_ws);   
                ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = ws_fig.colorbar(bias_ws,ax=ws_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                ws_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); ws_ax.set_yticks(HGTtick_hgh);ws_ax.set_yticklabels(HGTtick_hgh_label);
                ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                ws_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                ws_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = ws_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################                     
                                     
                """ ###############################################################
                Define and Plot the Model Wind Direction Biases
                ############################################################### """
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_wd = wd_ax.imshow(MDLwd*np.nan,extent=MDLgrid,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_wd);   
                wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = wd_fig.colorbar(bias_wd,ax=wd_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                ###################################################################
                wd_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); wd_ax.set_yticks(HGTtick_hgh);wd_ax.set_yticklabels(HGTtick_hgh_label);
                wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                #######################################################################
                # set figure y-label
                wd_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                wd_ax.set_xticks(Xtick); wd_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    wd_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                wd_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                wd_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                #######################################################################
                # Denote that Measurement Data was Not Available
                wd_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = wd_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                snr_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                #wd_ax.annotate('NOAA ESRL',xy = plt_param.IDpad_sptl,xycoords = 'figure fraction',\
                #                verticalalignment = 'bottom',horizontalalignment = 'right',\
                #                fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                wd_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                wd_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl, xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = snr_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_SNRhgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                snr_fig.savefig(FILout,dpi=300); plt.close('all');    
                
                
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
    
                    """ Extended Wind Speed Plots """
                    ws_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Wind Speed Time Histories
                    ##########################################
                    ws_fig.suptitle(str(frequency)+' MHz RWP WS',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    ws_rwp = rwp_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=WSmin_low_ext,vmax=WSmax_low_ext,\
                                         origin='lower',aspect='auto',cmap=cmap_ws);   
                    rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = ws_fig.colorbar(ws_rwp,ax=rwp_ax,ticks = WStick_low_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
                    rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    #######################################################################
                    # set figure y-label
                    rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    rwp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert stn name/alt as title
                    T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Wind Speed Time Histories
                    ##########################################
                    ws_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLws_ext,\
                                               vmin = WSmin_low_ext,vmax = WSmax_low_ext,shading='gouraud',cmap=cmap_ws)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_low_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###################################################################
                    # overlay barbs
                    mdl_ax.barbs(MDLtme_ext_mesh[BARBmdl_ext_Loi],MDLhgt_ext_mesh[BARBmdl_ext_Loi],MDLu_ext[BARBmdl_ext_Loi],\
                                MDLv_ext[BARBmdl_ext_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin_low,HGTmax_low); mdl_ax.set_yticks(HGTtick_low);mdl_ax.set_yticklabels(HGTtick_low_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' WS'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Speed Biases
                    ########################################################### """
                    bias_ws = bias_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin_low,HGTmax_low); bias_ax.set_yticks(HGTtick_low);bias_ax.set_yticklabels(HGTtick_low_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    #######################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                                    #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WSlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    ws_fig.savefig(FILout,dpi=300); plt.close('all');
      
                    """ Extended Wind Direction Plots """
                    wd_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Wind Direction Time Histories
                    ##############################################
                    wd_fig.suptitle(str(frequency)+' MHz RWP WD',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                                    fontsize=plt_param.TITLEfs);
                    ###############################################################
                    wd_rwp = rwp_ax.imshow(MDLwd_ext*np.nan,extent=MDLgrid_ext,vmin=WDmin_low_ext,\
                                         vmax=WDmax_low_ext,origin='lower',aspect='auto',cmap=cmap_wd_low_ext);    
                    rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = wd_fig.colorbar(wd_rwp,ax=rwp_ax,ticks = WDtick_low_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
                    rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    #######################################################################
                    # set figure y-label
                    rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    rwp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert stn name/alt as title
                    T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Wind Direction Time Histories
                    ##############################################
                    wd_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLwd_ext,vmin = WDmin_low_ext,\
                                        vmax = WDmax_low_ext,shading='gouraud',cmap=cmap_wd_low_ext);    
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_low_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###################################################################
                    # overlay barbs
                    mdl_ax.barbs(MDLtme_ext_mesh[BARBmdl_ext_Loi],MDLhgt_ext_mesh[BARBmdl_ext_Loi],MDLu_ext[BARBmdl_ext_Loi],\
                                MDLv_ext[BARBmdl_ext_Loi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin_low,HGTmax_low); mdl_ax.set_yticks(HGTtick_low);mdl_ax.set_yticklabels(HGTtick_low_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' WD'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Direction Biases
                    ########################################################### """
                    bias_wd = bias_ax.imshow(MDLwd_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin_low,HGTmax_low); bias_ax.set_yticks(HGTtick_low);bias_ax.set_yticklabels(HGTtick_low_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    #######################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                                    #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = wd_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WDlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    wd_fig.savefig(FILout,dpi=300); plt.close('all');
                                         
                    """ Extended Signal-to-Noise Ratio Plots """
                    #snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                    snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[3,1,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Signal-to-Noise Time Histories
                    ###############################################
                    snr_fig.suptitle(str(frequency)+' MHz RWP Range-Corrected SNR',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    snr_rwp = rwp_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=-40,vmax=10,\
                                         origin='lower',aspect='auto',cmap=cmap_snr);   
                    rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = snr_fig.colorbar(snr_rwp,ax=rwp_ax,ticks = np.arange(-40,10.1,10),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect3)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('SNR, dB', fontsize=plt_param.LABELfs);
                    ###############################################################
                    rwp_ax.set_ylim(HGTmin_low,HGTmax_low); rwp_ax.set_yticks(HGTtick_low);rwp_ax.set_yticklabels(HGTtick_low_label);
                    rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    #######################################################################
                    # set figure y-label
                    rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    rwp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert stn name/alt as title
                    T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################          
                                         
                    """ ###########################################################
                    Plot the Model Wind Speed Biases
                    ########################################################### """
                    bias_ws = ws_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(bias_ws,ax=ws_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    ws_ax.set_ylim(HGTmin_low,HGTmax_low); ws_ax.set_yticks(HGTtick_low);ws_ax.set_yticklabels(HGTtick_low_label);
                    ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    ws_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    ws_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    B = ws_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################                 
                                         
                    """ ###########################################################
                    Define and Plot the Model Wind Direction Biases
                    ########################################################### """
                    bias_wd = wd_ax.imshow(MDLwd_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(bias_wd,ax=wd_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    wd_ax.set_ylim(HGTmin_low,HGTmax_low); wd_ax.set_yticks(HGTtick_low);wd_ax.set_yticklabels(HGTtick_low_label);
                    wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    #######################################################################
                    # set figure y-label
                    wd_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    wd_ax.set_xticks(Xtick_ext); wd_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        wd_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    wd_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    wd_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    wd_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_low,HGTmax_low]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    B = wd_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    snr_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###############################################################
                    wd_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                                    #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                    wd_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = snr_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_SNRlow_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    snr_fig.savefig(FILout,dpi=300); plt.close('all');    
                    
                    ###############################################################
                    ###############################################################
                    ###############################################################
                    
                    """ Extended Wind Speed Plots """
                    ws_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Wind Speed Time Histories
                    ##########################################
                    ws_fig.suptitle(str(frequency)+' MHz RWP WS',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    ws_rwp = rwp_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=WSmin_hgh_ext,vmax=WSmax_hgh_ext,\
                                         origin='lower',aspect='auto',cmap=cmap_ws);   
                    rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = ws_fig.colorbar(ws_rwp,ax=rwp_ax,ticks = WStick_hgh_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
                    rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    #######################################################################
                    # set figure y-label
                    rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    rwp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert stn name/alt as title
                    T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Wind Speed Time Histories
                    ##########################################
                    ws_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLws_ext,\
                                               vmin = WSmin_hgh_ext,vmax = WSmax_hgh_ext,shading='gouraud',cmap=cmap_ws)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(ws_mdl,ax=mdl_ax,ticks = WStick_hgh_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###################################################################
                    # overlay barbs
                    mdl_ax.barbs(MDLtme_ext_mesh[BARBmdl_ext_Hoi],MDLhgt_ext_mesh[BARBmdl_ext_Hoi],MDLu_ext[BARBmdl_ext_Hoi],\
                                MDLv_ext[BARBmdl_ext_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); mdl_ax.set_yticks(HGTtick_hgh);mdl_ax.set_yticklabels(HGTtick_hgh_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' WS'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Speed Biases
                    ########################################################### """
                    bias_ws = bias_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(bias_ws,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); bias_ax.set_yticks(HGTtick_hgh);bias_ax.set_yticklabels(HGTtick_hgh_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    #######################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                    B = bias_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    ws_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###############################################################
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
                    #######################################################################
                    # add a logo
                    imax = ws_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WShgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    ws_fig.savefig(FILout,dpi=300); plt.close('all');
      
                    """ Extended Wind Direction Plots """
                    wd_fig, (rwp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Wind Direction Time Histories
                    ##############################################
                    wd_fig.suptitle(str(frequency)+' MHz RWP WD',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                                    fontsize=plt_param.TITLEfs);
                    ###############################################################
                    wd_rwp = rwp_ax.imshow(MDLwd_ext*np.nan,extent=MDLgrid_ext,vmin=WDmin_hgh_ext,\
                                         vmax=WDmax_hgh_ext,origin='lower',aspect='auto',cmap=cmap_wd_hgh_ext);   
                    rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = wd_fig.colorbar(wd_rwp,ax=rwp_ax,ticks = WDtick_hgh_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
                    rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    #######################################################################
                    # set figure y-label
                    rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    rwp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert stn name/alt as title
                    T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Wind Direction Time Histories
                    ##############################################
                    wd_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLwd_ext,vmin = WDmin_hgh_ext,\
                                        vmax = WDmax_hgh_ext,shading='gouraud',cmap=cmap_wd_hgh_ext);  
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(wd_mdl,ax=mdl_ax,ticks = WDtick_hgh_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Direction, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs, labelpad=10);
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); mdl_ax.set_yticks(HGTtick_hgh);mdl_ax.set_yticklabels(HGTtick_hgh_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # overlay barbs
                    mdl_ax.barbs(MDLtme_ext_mesh[BARBmdl_ext_Hoi],MDLhgt_ext_mesh[BARBmdl_ext_Hoi],MDLu_ext[BARBmdl_ext_Hoi],\
                                MDLv_ext[BARBmdl_ext_Hoi],linewidth = 0.35, length = 4.5,barbcolor='black')
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' WD'+MDLstr,{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Direction Biases
                    ########################################################### """
                    bias_wd = bias_ax.imshow(MDLwd_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(bias_wd,ax=bias_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Dir, $\mathrm{deg}^{\circ}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); bias_ax.set_yticks(HGTtick_hgh);bias_ax.set_yticklabels(HGTtick_hgh_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
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
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    B = bias_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    wd_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###############################################################
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
                    #######################################################################
                    # add a logo
                    imax = wd_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WDhgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    wd_fig.savefig(FILout,dpi=300); plt.close('all');
                                         
                    """ Extended Signal-to-Noise Ratio Plots """
                    #snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                    snr_fig, (rwp_ax,ws_ax,wd_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[3,1,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Lidar Signal-to-Noise Time Histories
                    ###############################################
                    snr_fig.suptitle(str(frequency)+' MHz RWP Range-Corrected SNR',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    snr_rwp = rwp_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=-40,vmax=10,\
                                         origin='lower',aspect='auto',cmap=cmap_snr);   
                    rwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = snr_fig.colorbar(snr_rwp,ax=rwp_ax,ticks = np.arange(-40,10.1,10),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect3)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('SNR, dB', fontsize=plt_param.LABELfs);
                    ###############################################################
                    rwp_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); rwp_ax.set_yticks(HGTtick_hgh);rwp_ax.set_yticklabels(HGTtick_hgh_label);
                    rwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    # set figure y-label
                    rwp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    rwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    rwp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert stn name/alt as title
                    T = rwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################          
                                         
                    """ ###########################################################
                    Plot the Model Wind Speed Biases
                    ########################################################### """
                    bias_ws = ws_ax.imshow(MDLws_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = ws_fig.colorbar(bias_ws,ax=ws_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Speed, m $\mathrm{s}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    ws_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); ws_ax.set_yticks(HGTtick_hgh);ws_ax.set_yticklabels(HGTtick_hgh_label);
                    ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    ws_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    ws_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    B = ws_ax.set_title('WS Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################                 
                                         
                    """ ###########################################################
                    Define and Plot the Model Wind Direction Biases
                    ########################################################### """
                    bias_wd = wd_ax.imshow(MDLwd_ext*np.nan,extent=MDLgrid_ext,vmin=-3,vmax=3,\
                                     origin='lower',aspect='auto',cmap=cmap_bias);   
                    wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = wd_fig.colorbar(bias_wd,ax=wd_ax,ticks = np.arange(-3,3.1,1),pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Wind Dir, ($\mathrm{deg}^{\circ}$)', fontsize=plt_param.LABELfs);
                    ###############################################################
                    wd_ax.set_ylim(HGTmin_hgh,HGTmax_hgh); wd_ax.set_yticks(HGTtick_hgh);wd_ax.set_yticklabels(HGTtick_hgh_label);
                    wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    wd_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    wd_ax.set_xticks(Xtick_ext); wd_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        wd_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    wd_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    wd_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    wd_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin_hgh,HGTmax_hgh]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    B = wd_ax.set_title('WD Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    snr_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###############################################################
                    wd_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                                    #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
                    wd_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = snr_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_SNRhgh_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    snr_fig.savefig(FILout,dpi=300); plt.close('all');    
                    ###############################################################
                    # Delete Extended Information to Preclude Incorrect Plotting
                    del MDLtime_ext
