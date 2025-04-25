import os,sys
#sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
from psldateutils import *
import matplotlib.dates as mdates
import os
import glob
import pandas
from itertools import groupby
import numpy as np
import xarray as xr
from datetime import datetime, timedelta, timezone
from pslobs_vars_dict import *
from psldateutils import *
from pslutils import *
from pslobs_vars_dict import *


def read_psl_ascii_rwp_data(stnid,datestring,frequency,datatype='Wind',timeavg='hourly',psldata_path='/psd2data/obs/realtime',stationDict=None):

    yj=ymd_to_yj(datestring)

    rpath='Radar'+str(frequency)

    fdir='Ww'+datatype
    if not timeavg=='hourly':
        fdir=fdir+'SubHourly'

    inpath=os.path.join(psldata_path,rpath,fdir,stnid.lower(),yj[:4],yj[4:],'*')

    listdir = sorted(glob.glob(inpath))
    obsvars = variables['rwp']

    if len(listdir):
        #%%
        #read in the ascii files 
        #the idea is to read in the different times as blocks using groupby
        #IN0 is a dict with the high-resolution mode data, IN1 is a dict with the low-resolution mode data
        #the script reads in the wind data or the temperature data from the RASS, defined with fdir (depending on which variable I want to plot)
        pulsewidths=[]
        numbeams=[]
        numranges=[]
        # mode 0 vars
        IN0={}
        IN0['time'] = []
        IN0['pulsewidth']=[]
        IN0['numrange']=[]
        # mode 1 vars
        IN1={}
        IN1['time']=[]
        IN1['pulsewidth']=[]
        IN1['numrange']=[]
        if 'Wind' in fdir:
            filevars=['height','wind_speed','wind_direction','qc','rv1','rv2','rv3','cnt1','cnt2','cnt3','snr1','snr2','snr3','qc1','qc2','qc3']
            outvars=['base_time', 'time_offset', 'height', 'wind_speed', 'wind_direction', 'wind_u', 'wind_v', 'wind_quality', 'radial_velocity', 'radial_velocity_quality', 'radial_velocity_navg', 'radial_velocity_snr' ]
            for fv in filevars:
                IN0[fv]=[]
                IN1[fv]=[]
        
        
            count0=0
            count1=0
            for l in range(0,len(listdir)):
                with open(listdir[l]) as f:
                   data = []
                   grps = groupby(f,key=lambda x:(x.strip().startswith(stnid.lower()) or x.strip().startswith(stnid.upper())))
        
                   for k,v in grps:
                       if k:
                            t=next(grps)[1]
                            tt=next(t)
                            tt=next(t) #lat, lon, alt
                            latlonalt=tt.split()
                            if not stationDict:
                              stationDict=dict(stn_id = stnid.upper(), lat = latlonalt[0], lon = latlonalt[1], alt=latlonalt[2])
                            tt=next(t) #date
                            dateinfo=tt.split()
                            #print(dateinfo)
                            utcoffset=timezone(timedelta(hours=int(dateinfo[6])))
                            dt=mdates.num2date(mdates.datestr2num(datestring[0:2]+dateinfo[0]+dateinfo[1]+dateinfo[2]+dateinfo[3]+dateinfo[4]+dateinfo[5]),tz=utcoffset)
                            #print(dt)
                            tt=next(t)  # avgtime, numbeams, numrange
                            numbeam=int(tt.split()[1])
                            if (numbeam != 3):
                               print("Number of beams != 3")
                               return None
                            numrange=int(tt.split()[2])
                            #print('num range',numrange)
                            tt=next(t)
                            tt=next(t) #pulse width
                            pulsewidth=int(tt.split()[4])
                            if (not pulsewidth in pulsewidths):
                               pulsewidths.append(pulsewidth)
                            tt=next(t)
                            tt=next(t)
                            tt=next(t)
                            # print(pulsewidth)
                            #if pulsewidth == pulse_width[0]:
                            if pulsewidth == pulsewidths[0]:
                                if count0==0:
                                    numrange0=numrange
                                    numranges.append(numrange0)
                                    numbeams.append(numbeam)
                                    for fv in filevars:
                                        IN0[fv]=[[] for j in range(numrange)]
                                        count0=1
                                if (numrange != numrange0):
                                   print('different mode 0 range',numrange0,numrange)
                                   return None
                                if dt in IN0['time']:
                                  print(dt,'already in list for pw0, replacing...')
                                  for i in range(0,numrange):
                                    tt=next(t)
                                    vals=tt.split()
                                    for i_fv in range(0,len(filevars)):
                                        IN0[filevars[i_fv]][i-numrange][0]=float(vals[i_fv])
                                else:
                                  #IN0['pulsewidth'].append(pulsewidth)
                                  IN0['time'].append(dt)
                                  for i in range(0,numrange):
                                    tt=next(t)
                                    vals=tt.split()
                                    for i_fv in range(0,len(filevars)):
                                        IN0[filevars[i_fv]][i].append(float(vals[i_fv]))
                            #elif pulsewidth == pulse_width[1]:
                            elif pulsewidth == pulsewidths[1]:
                                if count1==0:
                                    numrange1=numrange
                                    numranges.append(numrange1)
                                    numbeams.append(numbeam)
                                    for fv in filevars:
                                        IN1[fv]=[[] for j in range(numrange)]
                                        count1=1
                                if (numrange != numrange1):
                                   print('different mode 1 range',numrange1,numrange)
                                   return None
                                if dt in IN1['time']:
                                  print(dt,'already in list for pw1, replacing....')
                                  for i in range(0,numrange):
                                    tt=next(t)
                                    vals=tt.split()
                                    for i_fv in range(0,len(filevars)):
                                        IN1[filevars[i_fv]][i-numrange][0]=float(vals[i_fv])
                                else:
                                  #IN1['pulsewidth'].append(pulsewidth)
                                  IN1['time'].append(dt)
                                  for i in range(0,numrange):
                                      tt=next(t)
                                      vals=tt.split()
                                      for i_fv in range(0,len(filevars)):
                                          IN1[filevars[i_fv]][i].append(float(vals[i_fv]))
                            else: 
                                print('New pulsewidth',pulsewidth)

            numtimes=len(IN0['time'])
            IN0['time']=np.asarray(IN0['time'])
            IN1['time']=np.asarray(IN1['time'])
            #print('IN0 times = ',len(IN0['time']))
            #print('IN1 times = ',len(IN1['time']))
            #create base_time and time_offset from datetime
            #bt=datetime.strptime(datestring,'%Y%m%d').replace(tzinfo=timezone.utc)
            bt=mdates.num2date(mdates.datestr2num(datestring),tz=timezone.utc)
            IN0['base_time']=datetime.timestamp(bt)
            IN1['base_time']=datetime.timestamp(bt)
            #bt=bt.replace(tzinfo=None)
            IN0['time_offset']=[(dt-bt).total_seconds() for dt in IN0['time']]
            IN1['time_offset']=[(dt-bt).total_seconds() for dt in IN1['time']]

            for filevar in filevars:
                IN0[filevar]=np.transpose(np.asarray(IN0[filevar]))
                IN1[filevar]=np.transpose(np.asarray(IN1[filevar]))
                if ('wind_speed' in filevar) or ('wind_direction' in filevar):
                    IN0[filevar][IN0[filevar]==999999]=np.nan
                    IN1[filevar][IN1[filevar]==999999]=np.nan
            IN0['u'],IN0['v']=met2cart(IN0['wind_direction'],IN0['wind_speed'])
            IN1['u'],IN1['v']=met2cart(IN1['wind_direction'],IN1['wind_speed'])
            
        elif 'Temp' in fdir:
            filevars=['height','tv_uc','tv','w','qc_tv_uc','qc_tv','qc_w','cnt_tv_uc','cnt_tv','cnt_w','snr_tv_uc','snr_tv','snr_w']
            #outvars=['base_time','time_offset'] + filevars
            outvars=['base_time','time_offset'] 
            for fv in filevars:
                IN0[fv]=[]
        
            count0=0
            for l in range(0,len(listdir)):
                with open(listdir[l]) as f:
                   data = []
                   grps = groupby(f,key=lambda x:(x.strip().startswith(stnid.lower()) or x.strip().startswith(stnid.upper())))
        
                   for k,v in grps:
                       if k:
                            t=next(grps)[1]
                            tt=next(t) # data type; revision number
                            tt=next(t) #lat, lon, alt
                            tt=next(t) #date
                            dateinfo=tt.split()
                            utcoffset=timezone(timedelta(hours=int(dateinfo[6])))
                            dt=mdates.num2date(mdates.datestr2num(datestring[0:2]+dateinfo[0]+dateinfo[1]+dateinfo[2]+dateinfo[3]+dateinfo[4]+dateinfo[5]),tz=utcoffset)
                            tt=next(t)  # avg time, num beams, num range gates
                            numbeam=int(tt.split()[1])
                            numrange=int(tt.split()[2])
                            tt=next(t) # consensus
                            tt=next(t) #pulse width
                            pulsewidth=int(tt.split()[2])
                            if (not pulsewidth in pulsewidths):
                               pulsewidths.append(pulsewidth)
                            tt=next(t)
                            tt=next(t)
                            tt=next(t)
                            #if pulsewidth == pulse_width[0]:
                            if pulsewidth == pulsewidths[0]:
                                if count0==0:
                                    numrange0=numrange
                                    numbeams.append(numbeam)
                                    numranges.append(numrange0)
                                    for fv in filevars:
                                        IN0[fv]=[[] for j in range(numrange)]
                                        count0=1
                                if (numrange != numrange0):
                                   print('different range number ',numrange0,numrange)
                                   return None
                                if dt in IN0['time']:
                                  print(dt,'already in list for pw, replacing...')
                                  for i in range(0,numrange):
                                    tt=next(t)
                                    for i_fv in range(0,len(filevars)):
                                        IN0[filevars[i_fv]][i][i-numrange]=float(tt.split()[i_fv])
                                else:
                                  #IN0['pulsewidth'].append(pulsewidth)
                                  IN0['time'].append(dt)
                                  for i in range(0,numrange):
                                    tt=next(t)
                                    for i_fv in range(0,len(filevars)):
                                        IN0[filevars[i_fv]][i].append(float(tt.split()[i_fv]))
            IN0['time']=np.asarray(IN0['time'])
            #create base_time and time_offset from datetime
            bt=mdates.num2date(mdates.datestr2num(datestring),tz=timezone.utc)
            IN0['base_time']=datetime.timestamp(bt)
            IN0['time_offset']=[(dt-bt).total_seconds() for dt in IN0['time']]
            for filevar in filevars:
                IN0[filevar]=np.transpose(np.asarray(IN0[filevar]))
                if ('tv' in filevar) or ('w' in filevar):
                    IN0[filevar][IN0[filevar]==999999]=np.nan
        
        ## create the dataset
        OUT=xr.Dataset()
    
        ## create a time coordinate 
        timestamps=xr.DataArray(data=[datetime.timestamp(dt.replace(tzinfo=timezone.utc)) for dt in IN0['time']],dims='time')
        OUT['time']=timestamps
        OUT['time'].attrs= { 'long_name' : 'time',
                             'units' : 'seconds since 1970-01-01 00:00:00',
                             'calendar': 'standard',
                             'actual_range': [timestamps[0],timestamps[-1]]
                           }
        ## set time to be a double type
        OUT['time'].encoding.update(dict(dtype='f8'))
    
        ##  add in the station info to the dataset
        add_station_info(stationDict, OUT)
    
        modenums=[]
        nummodes=len(pulsewidths)
        for i in range(nummodes):
          modenums.append(i)
        OUT['mode_number'] = xr.DataArray(data=np.array(modenums),dims='mode')
        OUT['mode_number'].attrs= { 'long_name' : 'Mode number',
                             'units' : '1',
                             'flag_values' : np.array([0,1]),
                             'flag_meanings' : 'wind_low_power_mode wind_high_power_mode'
                           }
        OUT['num_range_gates'] = xr.DataArray(data=np.array(numranges),dims='mode')
        OUT['num_range_gates'].attrs= { 'long_name' : 'Number of range gates sampled per mode',
                           'units' : 1,
                           }
        OUT['num_beams'] = xr.DataArray(data=np.array(numbeams),dims='mode')
        OUT['num_beams'].attrs= { 'long_name' : 'Number of beam positions',
                           'units' : 1,
                           }
        OUT['pulse_width'] = xr.DataArray(data=np.array(pulsewidths),dims='mode')
        OUT['pulse_width'].attrs= { 'long_name' : 'Pulse width',
                             'units' : 'ns',
                             'actual_range': np.array([np.min(pulsewidths), np.max(pulsewidths)])
                           }

        # Create the height var - the levels is the max number of levels for all modes
        level=np.max(numranges)
        heightdata=np.empty((nummodes,level), dtype=np.float32);
        heightdata.fill(np.NAN)
        heightdata[0,:numranges[0]]=IN0['height'][0,:numranges[0]]
        if (nummodes>1):
           heightdata[1,:numranges[1]]=IN1['height'][0,:numranges[1]]
        #print(heightdata)
        d_=xr.DataArray(data=heightdata,dims=['mode','level'])
        OUT['height']=d_

        # merge the modes into a single var
        if (datatype=='Wind'):
            # wind vars
            windspeed=np.empty((numtimes,nummodes,level), dtype=np.float32)
            winddir=np.empty((numtimes,nummodes,level), dtype=np.float32)
            windqc=np.empty((numtimes,nummodes,level), dtype=np.int32)
            u=np.empty((numtimes,nummodes,level), dtype=np.float32)
            v=np.empty((numtimes,nummodes,level), dtype=np.float32)
            windspeed.fill(np.NAN)
            winddir.fill(np.NAN)
            windqc.fill(9)
            u.fill(np.NAN)
            v.fill(np.NAN)
            for time in range(numtimes):
                for mode in range(nummodes):
                    indict=IN0
                    if (mode==1):
                       indict=IN1
                    windspeed[time,mode,0:+numranges[mode]]=indict['wind_speed'][time]
                    winddir[time,mode,0:+numranges[mode]]=indict['wind_direction'][time]
                    windqc[time,mode,0:+numranges[mode]]=indict['qc'][time]
                    u[time,mode,0:+numranges[mode]]=indict['u'][time]
                    v[time,mode,0:+numranges[mode]]=indict['v'][time]
            OUT['wind_speed']=xr.DataArray(data=windspeed,dims=['time','mode','level'])
            OUT['wind_direction']=xr.DataArray(data=winddir,dims=['time','mode','level'])
            OUT['wind_quality']=xr.DataArray(data=windqc,dims=['time','mode','level'])
            OUT['wind_u']=xr.DataArray(data=u,dims=['time','mode','level'])
            OUT['wind_v']=xr.DataArray(data=v,dims=['time','mode','level'])
            # radial velocity
            rv=np.empty((numtimes,nummodes,3,level), dtype=np.float32)
            rv.fill(np.NAN)
            snr=np.empty((numtimes,nummodes,3,level), dtype=np.float32)
            snr.fill(np.NAN)
            cnt=np.zeros((numtimes,nummodes,3,level), dtype=np.float32)
            qc=np.zeros((numtimes,nummodes,3,level), dtype=np.float32)
            for time in range(numtimes):
                for mode in range(nummodes):
                    indict=IN0
                    if (mode==1):
                       indict=IN1
                    for beam in range(3):
                       if beam == 0: 
                           beamvar=indict['rv1']
                           snrvar=indict['snr1']
                           cntvar=indict['cnt1']
                           qcvar=indict['qc1']
                       elif beam == 1: 
                           beamvar=indict['rv2']
                           snrvar=indict['snr2']
                           cntvar=indict['cnt2']
                           qcvar=indict['qc2']
                       else:
                           beamvar=indict['rv3']
                           snrvar=indict['snr3']
                           cntvar=indict['cnt3']
                           qcvar=indict['qc3']
                       rv[time,mode,beam,0:+numranges[mode]]=beamvar[time]
                       snr[time,mode,beam,0:+numranges[mode]]=snrvar[time]
                       snr[snr==999999]=np.nan
                       cnt[time,mode,beam,0:+numranges[mode]]=cntvar[time]
                       qc[time,mode,beam,0:+numranges[mode]]=qcvar[time]
            OUT['radial_velocity']=xr.DataArray(data=rv,dims=['time','mode','beam','level'])
            OUT['radial_velocity_quality']=xr.DataArray(data=qc,dims=['time','mode','beam','level'])
            OUT['radial_velocity_navg']=xr.DataArray(data=cnt,dims=['time','mode','beam','level'])
            OUT['radial_velocity_snr']=xr.DataArray(data=snr,dims=['time','mode','beam','level'])
            

        for var in ['base_time','time_offset']:
            #print('writing out',var)
            if var=='base_time':
                d_=xr.DataArray(data=IN0[var])
            elif var=='time_offset':
                d_=xr.DataArray(data=IN0[var],dims='time')
            OUT[var]=d_
            OUT[var].encoding.update(dict(dtype=np.float32,zlib=True,complevel=2))   #  set these data variables to be float, compressed

        for var in outvars:
    
            ## add in the variable metadata
            min=np.min(OUT[var])     
            max=np.max(OUT[var])     

            obsvar=next((obsvar for obsvar in obsvars if obsvar['varname'] == var), None)
            if (obsvar):
                 varattrs=obsvar.get('attrs')
                 if '_FillValue' in varattrs.keys():
                      OUT[var].encoding.update(dict(_FillValue=varattrs['_FillValue']))   #  set the fill falue in encoding too 
                      varattrs.pop('_FillValue')
                 OUT[var].attrs=varattrs
            OUT[var].attrs['actual_range']=np.array([min, max], OUT[var].dtype)
    
        nowutc = datetime.utcnow()
        OUT.attrs= { 'Conventions' : 'CF-6.0',
                     'featureType' : 'timeSeriesProfile',
                     'creation_date' : nowutc.strftime('%Y-%m-%dT%H:%M:%SZ'),
                     'history': 'created '+nowutc.strftime('%Y/%m')+' by PSL Data & Web Group (non-packed netCDF4-Classic)\npsl.data@noaa.gov'
                   }
        return OUT
    else:
        print('no rwp'+str(frequency)+' data found at site '+stnid.upper()+ ' on '+datestring)
        return None

"""    
    #writes desired variable, time (as datenum), and height to dict OUT for plotting
    if 'snr' in var:
        plotinfo['unit'+var] = 'SNR (dB)'
    elif 'ffhm' in var:
        OUT[var]=IN0['wind_speed']
        plotinfo['unit'+var] = 'Wind speed (m s$^{-1}$)'
        #plotinfo['clim'+var] = [-3,3]
    elif 'fflm' in var:
        OUT[var]=IN1['wind_speed']
        plotinfo['unit'+var] = 'Wind speed (m s$^{-1}$)'
        #plotinfo['clim'+var] = [-3,3]
    elif 'windhm' in var:
        OUT['u'+var]=IN0['u']
        OUT['v'+var]=IN0['v']
    elif 'tvhm' in var:
        OUT[var]=IN0['tv']
        plotinfo['unit'+var] = 'Virtual temperature (deg C)'
    
    if 'hm' in var:
        OUT['height'+var] = IN0['height'][0,:]*1000
        OUT['datetime'+var] = IN0['time']
    elif 'lm' in var:
        OUT['height'+var] = IN1['height'][0,:]*1000
        OUT['datetime'+var] = IN1['time']
"""
