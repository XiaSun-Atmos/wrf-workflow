import os,sys
#sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
import numpy as np
import xarray as xr
import glob
import os
from datetime import datetime, timedelta, timezone
from pslobs_vars_dict import *
from psldateutils import *
from pslutils import *

def read_pslsfcmet_from_ascii(siteid,datestring,filevars=['pressure','temperature','relative_humidity','wind_speed','vecwind_speed',\
                'wind_direction','stdwind_direction','battery_voltage','solar_radiation','net_radiation','precipitation','maxwind_speed'],\
                outvars=None, stationDict=None,location='/psd2data/obs/realtime/CsiDatalogger/SurfaceMet/'):
    """ Read the PSL Campbell Scientific Data Logger ASCII file for a day from the realtime PSL archive
        and return an xarray Dataset that can be used to create a CF compliant netcdf file
    
    Parameters
    ----------
    siteid: str
       the site id  (e.g, 'ctd')
    datestring: str
       date to choose in yyyymmdd format (e.g., 20230401)
    filevars:  list       
       list of the input file variable names for each of the data variables 
    outvars:  list       
       list of the variables to save the output file (default=all)
    stationDict:  dict
        dictionary of station metadata (optional, but needed for CF compliance).  
        dictionary keys:   stn_name, stn_id, lat, lon, alt
    location: str
        top level directory for ascii files.  dir structure is location/siteid/julian_day
    
    Returns
    -------    
    ds  
        xarray Dataset.  Can be written out using ds.to_netcdf(filename).  Supports zlib compression if output is netcdf4.  
        See:  https://docs.xarray.dev/en/stable/generated/xarray.Dataset.to_netcdf.html
    """
    outvars = filevars if outvars is None else outvars
    ## sfcmet variable attributes
    obsvars = variables['sfc_met']
    ## sfcmet variable attributes (PSL uses lowercase ids)
    siteid = siteid.lower()
    fpathin=os.path.join(location,siteid)

    ## find the dates
    yj=ymd_to_yj(datestring)
    year=yj[:4]
    julday=yj[4:]
    #print(julday,year)
    l_ = sorted(glob.glob(os.path.join(fpathin,\
                                            year,julday,siteid+'*'+year[2:4]+julday+'.*m')))

    if len(l_):
        #print('process met data on '+ datestring)
        header=['id','year','jday','hmin'] + filevars

        IN={}
        for field in header:
            IN[field]=[]

        for l in range(0,len(l_)):
            #print('read '+l_[l])
            with open(l_[l]) as f:
                flist=list(f)
                # remove duplicate lines
                res = [i for n, i in enumerate(flist) if i not in flist[:n]]
                flist=res

            ## split on delimiters
            fsplit=[x.split(',') for x in flist]
            ffloat=list(np.float_(fsplit))
            data=np.asarray(ffloat)
            if data.shape[1] != len(header):
                print('wrong number of file variables in',l_[l])
                continue
            for field in header:
                idx=header.index(field)
                temp_=data[:,idx]
                # for some sites 400hPa has to be added to pressure - HACK! - won't work if name is not pressure
                if ('pressure' in field):
                    pmax=np.max(temp_)
                    if (pmax > -99999. and pmax < 670):
                       temp_ = temp_+400
                temp_[temp_ == -99999.] = 99999.
                IN[field]=np.hstack((IN[field],temp_))
        #fill hour with leading 0
        h=[f"{int(x):04}" for x in IN['hmin']]
        #print(h)
        jd=[str(int(x)) for x in IN['jday']]
        #convert to datetime
        IN['datetime']=np.asarray([datetime(int(IN['year'][i]),1,1) + timedelta(days=float(jd[i])-1,\
            hours=float(h[i][0:2]),minutes=float(h[i][2:])) for i in range(0,len(IN['year']))])
        #print(IN['datetime'])
        IN['base_time']=datetime.timestamp(IN['datetime'][0])


        #create base_time and time_offset from datetime
        bt=datetime.strptime(datestring,'%Y%m%d').replace(tzinfo=timezone.utc)
        IN['base_time']=datetime.timestamp(bt)
        bt=bt.replace(tzinfo=None)
        IN['time_offset']=[(dt-bt).total_seconds() for dt in IN['datetime']]


        ## create the dataset
        OUT=xr.Dataset()

        ## create a time coordinate 
        timestamps=xr.DataArray(data=[datetime.timestamp(dt.replace(tzinfo=timezone.utc)) for dt in IN['datetime']],dims='time')
        #print(timestamps)
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

        vars=['base_time','time_offset'] + outvars

        for var in vars:
            if var=='base_time':
                d_=xr.DataArray(data=IN[var])
            else:
                d_=xr.DataArray(data=IN[var],dims='time')
            OUT[var]=d_
            OUT[var].encoding.update(dict(dtype=np.float32,zlib=True,complevel=2))   #  set data variables to be float, compressed

            ## add in the variable metadata
            obsvar=next((obsvar for obsvar in obsvars if obsvar['varname'] == var), None)
            if (obsvar):
                 OUT[var].attrs=obsvar.get('attrs')
            min=np.min(OUT[var])     
            max=np.max(OUT[var])     
            OUT[var].attrs['actual_range']=np.array([min, max], np.float32)

        nowutc = datetime.utcnow()
        OUT.attrs= { 'Conventions' : 'CF-6.0',
                     'featureType' : 'timeSeries',
                     'creation_date' : nowutc.strftime('%Y-%m-%dT%H:%M:%SZ'),
                     'history': 'created '+nowutc.strftime('%Y/%m')+' by PSL Data & Web Group (non-packed netCDF4-Classic)\npsl.data@noaa.gov'
                   }
        return OUT
    else:
        print('no met data found at site '+siteid+ ' on '+datestring)
        return None


####def add_station_info(stationInfo, ds):
####    """ Add the station data (name, id, lat, lon and alt) to the dataset.  
####        The dictionary keys are: (stn_name, stn_id, lat, lon, alt)
####    """
####    if (stationInfo):
####        stnname = stationInfo.get('stn_name')
####        if (stnname):
####            ds['station_name'] = stnname
####            ds['station_name'].attrs = { 'long_name' : 'station_name'}
####            ds['station_name'].encoding.update(dict(dtype='S1'))
####        stid = stationInfo.get('stn_id')
####        if (stid):
####            ds['station_id'] = stid
####            ds['station_id'].attrs = { 'long_name' : 'station_id'}
####            ds['station_id'].encoding.update(dict(dtype='S1'))
####        lat = stationInfo.get('lat')
####        if (lat):
####            ds['lat'] = float(lat)
####            ds['lat'].attrs = { 'long_name' : 'station latitude',
####                                'units'     : 'degrees_north',
####                                'standard_name' : 'latitude'
####                               }
####            ds['lat'].encoding.update(dict(dtype='f4'))
####        lon = stationInfo.get('lon')
####        if (lon):
####            ds['lon'] = float(lon)
####            ds['lon'].attrs = { 'long_name' : 'station longitude',
####                                'units'     : 'degrees_west',    ## PSL standard is degrees_W
####                                'standard_name' : 'latitude'
####                               }
####            ds['lon'].encoding.update(dict(dtype='f4'))
####        alt = stationInfo.get('alt')
####        if (alt):
####            ds['alt'] = float(alt)
####            ds['alt'].attrs = { 'long_name' : 'station altitude',
####                                'units'     : 'm',
####                                'standard_name' : 'height'
####                               }
####            ds['alt'].encoding.update(dict(dtype='f4'))
####    return


