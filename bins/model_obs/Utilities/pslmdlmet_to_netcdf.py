import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__),'../../Utilities/'))
import numpy as np
import xarray as xr
import glob
from datetime import datetime, timedelta, timezone
from mdlobs_vars_dict import *
from psldateutils import *
from pslutils import *

def read_pslmdlmet_from_ascii(siteid,datestring,model,filevars=['precipitation'],outvars=None,stationDict=None,outputDir=None,modelDict=None):
    """ Read the PSL Campbell Scientific Data Logger ASCII file for a day from the realtime PSL archive
        and return an xarray Dataset that can be used to create a CF compliant netcdf file
    
    Parameters
    ----------
    siteid: str
       the site id  (e.g, 'ctd')
    datestring: str
       date to choose in yyyymmdd format (e.g., 20230401)
    model: str
       HRRR_v4 or RAP_v5
    filevars:  list       
       list of the input file variable names for each of the data variables 
    outvars:  list       
       list of the variables to save the output file (default=all)
    stationDict:  dict
        dictionary of station metadata (optional, but needed for CF compliance).  
        dictionary keys:   stn_name, stn_id, lat, lon, alt
    outputDir:  str
        path to the output directory to write the files
    modelDict:  dict
        dictionary of model grid point metadata 
        dictionary keys:   mlat, mlon, malt
    
    
    """
    outvars = filevars if outvars is None else outvars
    outdir = '.' if outputDir is None else outputDir
    obsvars = variables['mdl_met']
    ## pslmet variable attributes (PSL uses lowercase ids)
    siteid = siteid.lower()
    mdlid = model.lower();
    if ('hrrr' in mdlid):
       location='/psd2data/obsarchive/HRRR/3kmCONUS/'
    elif ('rap' in mdlid):
       location='/psd2data/obs/realtime/RAP/13kmCONUS/'
    fpathin=os.path.join(location,siteid)
    hrrr_forecast_hours=19
    hrrr_forecast_hours_ext=49
    rap_forecast_hours=22
    rap_forecast_hours_ext=52

    ## find the dates
    yj=ymd_to_yj(datestring)
    year=yj[:4]
    julday=yj[4:]
    #print(julday,year)
    inpath=os.path.join(fpathin,year,julday,siteid+year[2:4]+julday+'??.met')
    #print(inpath)
    l_ = sorted(glob.glob(inpath))

    if len(l_):
        #print('process met data on '+ datestring)
        header=['id','year','jday','hmin'] + filevars

        for l in range(0,len(l_)):
        #for l in range(1):

            IN={}
            for field in header:
                IN[field]=[]

            #print('read '+l_[l])
            with open(l_[l]) as f:
                flist=list(f)
                # remove duplicate lines
                res = [i for n, i in enumerate(flist) if i not in flist[:n]]
                flist=res

            ## Extract the run time from the filename.
            fhour=l_[l][-6:-4]
            fhour4=fhour+'00'
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
                #  Set any of the instrument malfunction vals (-99999) to the missing value (99999)
                temp_[temp_ == -99999.] = 99999.
                IN[field]=np.hstack((IN[field],temp_))
            #fill hour with leading 0
            h=[f"{int(x):04}" for x in IN['hmin']]
            #print(h)
            jd=[str(int(x)) for x in IN['jday']]
            #convert to datetime
            IN['datetime']=np.asarray([datetime(int(IN['year'][i]),1,1) + timedelta(days=float(jd[i])-1,\
                hours=float(h[i][0:2]),minutes=float(h[i][2:])+1) for i in range(0,len(IN['year']))])
            #print(IN['datetime'])
            #IN['base_time']=datetime.timestamp(IN['datetime'][0])
    
    
            #create base_time and time_offset from datetime
            bt=datetime.strptime(datestring+fhour,'%Y%m%d%H').replace(tzinfo=timezone.utc)
            IN['base_time']=datetime.timestamp(bt)
            bt=bt.replace(tzinfo=None)
            #IN['datetime']=np.insert(IN['datetime'],[0],[bt]);
            #bt=bt.replace(tzinfo=None)
            IN['time_offset']=[(dt-bt).total_seconds() for dt in IN['datetime']]
            IN['forecast']=[int(x/3600) for x in IN['time_offset']]
            numfcsts_in_file=len(IN['forecast'])
            numfcsts=numfcsts_in_file
            if ('rap' in model):
               if int(fhour) in [3,9,15,21]:
                  numfcsts=rap_forecast_hours_ext
               else:
                  numfcsts=rap_forecast_hours
            elif ('hrrr' in model):
               if int(fhour) in [0,6,12,18]:
                  numfcsts=hrrr_forecast_hours_ext
               else:
                  numfcsts=hrrr_forecast_hours
            fhours=np.arange(numfcsts)
            times=[]
            for i in range(numfcsts):
               times.append(datetime(int(year),1,1) + timedelta(days=int(julday)-1,\
                 hours=int(fhour)) + timedelta(hours=int(fhours[i])))
            IN['datetime']=np.asarray(times)
            #print(IN['forecast'])
            #forecast_hours=len(IN['forecast'])

    
    
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
    
            ##  add in the model info to the dataset
            add_model_station_info(modelDict, OUT)

            vars=['base_time','forecast'] + outvars
    
            for var in vars:
                #print(var)
                
                if var=='base_time':
                    d_=xr.DataArray(data=IN[var])
                elif var=='forecast':
                    d_=xr.DataArray(data=fhours,dims='forecast_hours')
                else:
                    if (len(IN[var]) < numfcsts):
                      temp=np.empty(numfcsts,dtype=np.float32)
                      temp.fill(np.nan)
                      ##indices=np.where(np.in1d(fhours,IN['forecast']))[0]
                      indices=[i for i, x in enumerate(fhours) if x in IN['forecast']]
                      #print(indices)
                      for idx in range(len(IN[var])):
                          temp[indices[idx]]=IN[var][idx]
                      #for i in range(len(fhours)):
                      #   if np.any(fhours[i]==IN['forecast']):
                      #      #print(fhours[i])
                      #      id=np.squeeze(np.argwhere(fhours[i]==IN['forecast']))
                      #      #print('pos',id)
                      #      temp[i]=IN[var][id]
                      IN[var]=temp
                    d_=xr.DataArray(data=IN[var],dims='forecast_hours')
                OUT[var]=d_
                OUT[var].encoding.update(dict(dtype=np.float32,zlib=True,complevel=2))   #  set data variables to be float, compressed
    
                ## add in the variable metadata
                obsvar=next((obsvar for obsvar in obsvars if obsvar['varname'] == var), None)
                if (obsvar):
                     OUT[var].attrs=obsvar.get('attrs')
                min=np.min(OUT[var])     
                max=np.max(OUT[var])     
                if (obsvar):
                   varattrs=obsvar.get('attrs')
                   if '_FillValue' in varattrs.keys():
                        OUT[var].encoding.update(dict(_FillValue=varattrs['_FillValue']))   #  set the fill falue in encoding too 
                        varattrs.pop('_FillValue')
                   OUT[var].attrs=varattrs
                OUT[var].attrs['actual_range']=np.array([min, max], OUT[var].dtype)
    
            nowutc = datetime.utcnow()
            OUT.attrs= { 'Conventions' : 'CF-6.0',
                         'featureType' : 'timeSeries',
                         'creation_date' : nowutc.strftime('%Y-%m-%dT%H:%M:%SZ'),
                         'history': 'created '+nowutc.strftime('%Y/%m')+' by PSL Data & Web Group (non-packed netCDF4-Classic)\npsl.data@noaa.gov',
                         'source': 'generated from '+ l_[l]
                       }
            outfile='mdlmet_'+siteid+'.'+model+'.'+datestring+'.'+fhour4+'.nc'
            #print("writing to",outfile)
            OUT.to_netcdf(os.path.join(outdir,outfile), format="NETCDF4_CLASSIC", engine='netcdf4')
    else:
        print('no model data found at site '+siteid+ ' on '+datestring)


