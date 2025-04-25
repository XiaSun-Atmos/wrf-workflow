import os,sys
#sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
import numpy as np
import xarray as xr
import glob
from datetime import datetime, timedelta, timezone
from pslobs_vars_dict import *
from psldateutils import *
from pslutils import *

#MM:SS:mmm MM:SS:mmm    B1    B2    B3    B4    B5    B6    B7    B8    B9   B10   B11   B12   B13   B14   B15   B16   B17   B18   B19   B20   B21   B22   B23   B24   B25   B26   B27   B28   B29   B30   B31   B32 | Blackout  Good   Bad | NumParticle  Rate(mm/h)  Amount(mm) AmountSum(mm)   Z(dB) | NumError Dirty VeryDirty Damaged SignalAvg SignalStdDev | TempAvg(C) TempStdDev(C) VoltAvg(V) VoltStdDev(V) HeatCurrentAvg(A) HeatCurrentStdDev(A) |   NumRain  NumNoRain  NumAmbig  Type
#00:00:000-02:00:000     0     0    17    49    43    21    19    21    25    30    81    65    36    10     8     8     7     1     1     0     0     0     0     0     0     0     0     0     0     0     0     0          0    12     0           442       0.570        0.02         27.42  16.816          0     0         0       0   24559.4        3.602         -2.0         0.000      21.00         0.000              3.92                0.002           0        350        75     3


def read_psldisdrometerstats_from_ascii(siteid,datestring, outvars=None, stationDict=None,location='/psd2data/obs/realtime/DisdrometerParsivel/Stats/'):

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
    filevars=['B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13','B14','B15','B16','B17','B18','B19','B20','B21','B22','B23','B24','B25','B26','B27','B28','B29','B30','B31','B32','Blackout','Good','Bad','NumParticle','Rate','Amount','AmountSum','Z','NumError','Dirty','VeryDirty','Damaged','SignalAvg','SignalStdDev','TempAvg','TempStdDev','VoltAvg','VoltStdDev','HeatCurrentAvg','HeatCurrentStdDev','NumRain','NumNoRain','NumAmbig','Type']
    outvars = filevars if outvars is None else outvars
    ## sfcmet variable attributes
    obsvars = variables['disdrometer']
    ## sfcmet variable attributes (PSL uses lowercase ids)
    siteid = siteid.lower()
    fpathin=os.path.join(location,siteid)
    #print(fpathin)

    ## find the dates
    yj=ymd_to_yj(datestring)
    year=yj[:4]
    julday=yj[4:]
    #print(julday,year)
    files_ = sorted(glob.glob(os.path.join(fpathin,\
                                            year,julday,siteid+'*'+year[2:4]+julday+'*_stats.txt')))

    if len(files_):
        #print('process met data on '+ datestring)
        header=['daterange'] + filevars

        IN={}
        for field in header:
            IN[field]=[]
        IN['hr_start']=[]

        for l in range(0,len(files_)):
            #print('read '+files_[l])
            with open(files_[l]) as f:
                flist=list(f)
                # remove duplicate lines
                res = [i for n, i in enumerate(flist) if i not in flist[:n]]
                flist=res
            if (len(flist)<4):
                continue

            ## Get the date from the first line
            dateline=flist[0].split()
            datestr=dateline[6]
            #print(datestr)
            ## split on delimiters
            fsplit=[x.split() for x in flist[3:]]
            #ffloat=list(np.float_(fsplit))
            #data=np.asarray(ffloat)
            data=np.asarray(fsplit)
            if data.shape[1] != len(header):
                print('wrong number of file variables in',files_[l])
                continue
            for field in header:
                idx=header.index(field)
                temp_=data[:,idx]
                temp=np.asarray(temp_)
                IN[field]=np.hstack((IN[field],temp))
            for h in range(len(temp)):
                IN['hr_start'].append(datestr[5:])
                
        jd=datestr[2:5]
        jdays=float(jd)-1
        #print('jd',jd)
        fyear='20'+datestr[:2]
        hr_start=datestr[5:]
        #print(IN['daterange'])
        dateranges=list(IN['daterange'])
        startend=[x.split('-') for x in dateranges]
        IN['starttimes'] = []
        IN['endtimes'] = []
        for t in range(len(startend)):
           IN['starttimes'].append(startend[t][0])
           IN['endtimes'].append(startend[t][1])
        #print(IN['starttimes'])
        #print(IN['endtimes'])
        #convert to datetime
        IN['datetimestart']=np.asarray([datetime(int(year),1,1) + timedelta(days=jdays,\
            hours=float(IN['hr_start'][i]),minutes=float(IN['starttimes'][i][:2]),seconds=float(IN['starttimes'][i][3:5]),milliseconds=float(IN['starttimes'][i][6:])) for i in range(0,len(IN['starttimes']))])
        #print(IN['datetimestart'][:40])
        IN['datetimeend']=np.asarray([datetime(int(year),1,1) + timedelta(days=jdays,\
            hours=float(IN['hr_start'][i]),minutes=float(IN['endtimes'][i][:2]),seconds=float(IN['endtimes'][i][3:5]),milliseconds=float(IN['endtimes'][i][6:])) for i in range(0,len(IN['endtimes']))])
        ##  account for last time in a file being the time on the next hour
        for dt in range(1,len(IN['datetimeend'])):
           if (IN['datetimeend'][dt] < IN['datetimeend'][dt-1]):
              doy=jdays
              offdt=IN['datetimeend'][dt]
              offhour=offdt.hour+1
              if (offhour >= 24):
                 doy=doy+1
                 offhour=0
              IN['datetimeend'][dt]=datetime(offdt.year,1,1)+timedelta(days=doy,hours=offhour,minutes=offdt.minute,seconds=offdt.second,microseconds=offdt.microsecond)

        #print(IN['datetimeend'][:40])
        #print(IN['datetime'])
        #IN['base_time']=datetime.timestamp(IN['datetime'][0])
        #IN['base_time']=datetime.timestamp(datetime(int(year),int(datestring[4:6]),int(datestring[6:])))


        #create base_time and time_offset from datetime
        bt=datetime.strptime(datestring,'%Y%m%d').replace(tzinfo=timezone.utc)
        IN['base_time']=datetime.timestamp(bt)
        bt=bt.replace(tzinfo=None)
        IN['time_offset']=[(dt-bt).total_seconds() for dt in IN['datetimestart']]
        #print(IN['time_offset'])

        #calculate the observation duration
        IN['time_interval']=[dt.total_seconds() for dt in (IN['datetimeend']-IN['datetimestart'])]
####
####
####        ## create the dataset
        OUT=xr.Dataset()

        ## create a time coordinate 
        timestamps=xr.DataArray(data=[datetime.timestamp(dt.replace(tzinfo=timezone.utc)) for dt in IN['datetimestart']],dims='time')
        timestampsend=xr.DataArray(data=[datetime.timestamp(dt.replace(tzinfo=timezone.utc)) for dt in IN['datetimeend']],dims='time')
        #print(timestamps)
        OUT['time']=timestamps
        OUT['time'].attrs= { 'long_name' : 'time',
                             'units' : 'seconds since 1970-01-01 00:00:00',
                             'calendar': 'standard',
                             'actual_range': [timestamps[0],timestamps[-1]]
                           }
        OUT['time_end']=timestampsend
        OUT['time_end'].attrs= { 'long_name' : 'observation end time',
                             'units' : 'seconds since 1970-01-01 00:00:00',
                             'calendar': 'standard',
                             'actual_range': [timestampsend[0],timestampsend[-1]]
                           }
        ## set time to be a double type
        OUT['time'].encoding.update(dict(dtype='f8'))
        OUT['time_end'].encoding.update(dict(dtype='f8'))


        ##  add in the station info to the dataset
        add_station_info(stationDict, OUT)
        outvars=['Rate','Amount','AmountSum']

        vars=['base_time','time_offset', 'time_interval'] + outvars

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
        print('no disdrometer data found at site '+siteid+ ' on '+datestring)
        return None


