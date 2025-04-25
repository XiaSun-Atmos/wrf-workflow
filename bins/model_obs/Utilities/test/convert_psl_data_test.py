#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
###############################################################################
Created on Tue Oct 22 11:31:51 2019
@author: don.murray
###############################################################################
"""
import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__),'../'))
from datetime import date,timedelta
from pslsfcmet_to_netcdf import read_pslsfcmet_from_ascii
from psldisd_to_netcdf import read_psldisdrometerstats_from_ascii
from pslrwp_to_netcdf import read_psl_ascii_rwp_data


def DATAextr_auto(Soi,DSsrch,Stime,Etime,DIRout):
    ###########################################################################
    #Establish Output Directory    
    output_dir = os.path.join(DIRout,Soi+'out')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ###########################################################################
    ### Identify Date Range to Examine
    ###########################################################################
    # Identify the Period of Interest
    strt_date = Stime
    end_date = Etime
    #print(Stime,Etime)
    ###########################################################################
    # Define the Date Range of Interest
    anly_dates = []
    days = (date(int(end_date[0:4]),int(end_date[4:6]),int(end_date[6:8]))\
          -date(int(strt_date[0:4]),int(strt_date[4:6]),int(strt_date[6:8]))).days
    drnge = range(days+1)
    for d in drnge:
        dt = date(int(strt_date[0:4]),int(strt_date[4:6]),int(strt_date[6:8]))+timedelta(d)
        #######################################################################
        anly_dates.append(dt.strftime('%Y%m%d'))

    ###########################################################################
    # Process the raw text data to netCDF for the specified dates
    """ Perform Relevant PSL Measurement Data Storage
    ####################################################################### """
    # Identify Relevant Measurement Datastreams to Extract (Station Specific)
    for dtype in DSsrch:
        for d2prc in anly_dates:
            """ Give Update on Data Storage Output
            ################################################################### """
            print('Storing ' + dtype + ' for ' + Soi + ':: ' + d2prc)
            """ Initiate NetCDF File """
            if dtype == 'pslmet':
                DSfilname = os.path.join(output_dir,dtype + '.' + Soi + '.' + d2prc + '.nc')
                if (os.path.exists(DSfilname)):
                    os.remove(DSfilname)
                xrf=read_pslsfcmet_from_ascii(Soi.lower(),d2prc, stationDict=stninfo, filevars=fvars)
                if xrf != None:
                    xrf.to_netcdf(DSfilname, format="NETCDF4_CLASSIC", engine='netcdf4')
            elif 'rwp915' in dtype:
                datatype='Wind'
                freq = '915'
                timetypes = ['hourly','subhourly']
                for timeavg in timetypes:
                   #rwp915windhourly.CTD.20230101.nc
                   #sgp915rwpwindavgC1.a1.20210406.100000.nc
                   DSfilname = os.path.join(output_dir,'rwp' + freq + datatype.lower() + timeavg + '.' + Soi + '.' + d2prc + '.nc')
                   if (os.path.exists(DSfilname)):
                       os.remove(DSfilname)
                   xrf = read_psl_ascii_rwp_data(Soi.lower(),d2prc,freq,datatype=datatype,timeavg=timeavg,psldata_path='/psd2data/obs/realtime',stationDict=stninfo)
                   if xrf != None:
                       xrf.to_netcdf(DSfilname, format="NETCDF4_CLASSIC", engine='netcdf4')

            elif 'rwp449' in dtype:
                datatype='Wind'
                freq = '449'
                timetypes = ['hourly','subhourly']
                for timeavg in timetypes:
                    #rwp449windhourly.CTD.20230101.nc
                    DSfilname = os.path.join(output_dir,'rwp' + freq + datatype.lower() + timeavg + '.' + Soi + '.' + d2prc + '.nc')
                    if (os.path.exists(DSfilname)):
                        os.remove(DSfilname)
                    xrf = read_psl_ascii_rwp_data(Soi.lower(),d2prc,freq,datatype=datatype,timeavg=timeavg,psldata_path='/psd2data/obs/realtime',stationDict=stninfo)
                    if xrf != None:
                        xrf.to_netcdf(DSfilname, format="NETCDF4_CLASSIC", engine='netcdf4')

            elif (dtype == 'psldisd'):
                DSfilname = os.path.join(output_dir,dtype + '.' + Soi + '.' + d2prc + '.nc')
                xrf=read_psldisdrometerstats_from_ascii(Soi.lower(),d2prc,stationDict=stninfo)
                if xrf != None:
                    xrf.to_netcdf(DSfilname, format="NETCDF4_CLASSIC", engine='netcdf4')
            else:
                continue

####            have_new_data = False;
####            for ds_oi in range(0,len(Soi_dstream)):
####                # Define the Station/Datastream Measurement Files #
####                Sds_fils = [f for f in os.listdir(meas_dir) if f.startswith(Soi_dstream[ds_oi])]
####                # Identify File Associated with Date
####                Sds_foi = [f for f in Sds_fils if d2prc in f];
####                #print(ds_oi)
####                #print(Sds_foi)
####                any_new = check_for_new_data(DSfilname, Sds_foi, meas_dir)
####                #print('any new? '+str(any_new))
####                have_new_data = have_new_data or any_new
####            #print('have new data? '+str(have_new_data))
####    
####            if not have_new_data:
####               #print('No new data for '+d2prc)
####               continue
####            else:
####               print('Have new data for '+d2prc)


# Run with hard coded station information
Soi='CTD'
DSrch = ['pslmet','rwp915','rwp449','psldisd']
DIRout='.'
stninfo = {
    'stn_id'  : 'CTD',
    'stn_name' : 'Courtland, AL',
    'lat'      : 34.66,
    'lon'      : -87.35,
    'alt'      : 187,
}
fvars = ['pressure','temperature','relative_humidity','wind_speed',\
         'vecwind_speed','wind_direction','stdwind_direction',\
         'battery_voltage','solar_radiation','net_radiation',\
         'precipitation','maxwind_speed'];
Stime='20230401'
Etime='20230403'

DATAextr_auto(Soi,DSrch,Stime,Etime,DIRout)
