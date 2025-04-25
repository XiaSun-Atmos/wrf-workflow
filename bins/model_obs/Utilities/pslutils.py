import numpy as np
import os

def add_station_info(stationInfo, ds):
    """ Add the station data (name, id, lat, lon and alt) to the dataset.  
        The dictionary keys are: (stn_name, stn_id, lat, lon, alt)
    """
    if (stationInfo):
        stnname = stationInfo.get('stn_name')
        if (stnname):
            ds['station_name'] = stnname
            ds['station_name'].attrs = { 'long_name' : 'station_name'}
            ds['station_name'].encoding.update(dict(dtype='S1'))
        stid = stationInfo.get('stn_id')
        if (stid):
            ds['station_id'] = stid
            ds['station_id'].attrs = { 'long_name' : 'station_id'}
            ds['station_id'].encoding.update(dict(dtype='S1'))
        lat = stationInfo.get('lat')
        if (lat):
            ds['lat'] = float(lat)
            ds['lat'].attrs = { 'long_name' : 'station latitude',
                                'units'     : 'degrees_north',
                                'standard_name' : 'latitude'
                               }
            ds['lat'].encoding.update(dict(dtype='f4'))
        lon = stationInfo.get('lon')
        if (lon):
            ds['lon'] = float(lon)
            ds['lon'].attrs = { 'long_name' : 'station longitude',
                                'units'     : 'degrees_west',    ## PSL standard is degrees_W
                                'standard_name' : 'latitude'
                               }
            ds['lon'].encoding.update(dict(dtype='f4'))
        alt = stationInfo.get('alt')
        if (alt):
            ds['alt'] = float(alt)
            ds['alt'].attrs = { 'long_name' : 'station altitude',
                                'units'     : 'm',
                                'standard_name' : 'height'
                               }
            ds['alt'].encoding.update(dict(dtype='f4'))
    return

def add_model_station_info(modelInfo, ds):
    """ Add the model station grid point location data (lat, lon and alt) to the dataset.  
        The dictionary keys are: (mlat, mlon, malt)
    """
    if (modelInfo):
        mlat = modelInfo.get('mlat')
        if (mlat):
            ds['mlat'] = float(mlat)
            ds['mlat'].attrs = { 'long_name' : 'model grid point latitude',
                                'units'     : 'degrees_north',
                                'standard_name' : 'latitude'
                               }
            ds['mlat'].encoding.update(dict(dtype='f4'))
        mlon = modelInfo.get('mlon')
        if (mlon):
            ds['mlon'] = float(mlon)
            ds['mlon'].attrs = { 'long_name' : 'model grid point longitude',
                                'units'     : 'degrees_east',
                                'standard_name' : 'latitude'
                               }
            ds['mlon'].encoding.update(dict(dtype='f4'))
        malt = modelInfo.get('malt')
        if (malt):
            ds['malt'] = float(malt)
            ds['malt'].attrs = { 'long_name' : 'model grid point altitude',
                                'units'     : 'm',
                                'standard_name' : 'height'
                               }
            ds['malt'].encoding.update(dict(dtype='f4'))
    return

                
def met2cart(dd,ff):
	#transforms meteorlogical wind direction ddinto cartesian components u,v
	u=-ff*np.sin(dd*np.pi/180.)
	v=-ff*np.cos(dd*np.pi/180.)

	return u,v

def get_tropoe_info(project,TROPoetree,tropoe_dir,tropoe_sid,tropoe_pid,version):
    if 'perils' in project:
        suffix = '.cdf'
        if 'assist' in version:
          sid_path = os.path.join(TROPoetree,tropoe_dir,'retrieval_output');
          f_pfx=tropoe_pid.lower()+'tropoe'+tropoe_dir.lower()+'.v'
          dateindex=2
        else:
          sid_path = os.path.join(TROPoetree,tropoe_dir,'retrieval_output','zo');
          f_pfx=tropoe_pid.lower()+'tropoe'+tropoe_dir+'.c1.zo'
          dateindex=3
    elif 'wfip3' in project:
        suffix = '.nc'
        if 'assist' in version:
          sid_path = os.path.join(TROPoetree,tropoe_dir,'retrieval_output');
          f_pfx=tropoe_pid.lower()+'tropoe'+tropoe_dir+'.rt.2'
          dateindex=2
        else:
          sid_path = os.path.join(TROPoetree,tropoe_dir,'retrieval_output');
          f_pfx=tropoe_pid.lower()+'tropoe'+tropoe_dir+'.rt.2'
          dateindex=2
    elif 'splash' in project:
        suffix = '.cdf'
        if 'assist' in version:
          sid_path = os.path.join(TROPoetree,tropoe_dir+'_realtime','retrieval_output');
          f_pfx=tropoe_pid.lower()+'tropoe'+tropoe_dir.lower()+'.v'
          dateindex=2
        else:
          sid_path = os.path.join(TROPoetree,tropoe_dir,'retrieval_output');
          f_pfx=tropoe_sid.lower()+'tropoe_'+tropoe_dir.lower()+'.v'
          dateindex=2
    else:
        suffix = '.nc'
        sid_path = os.path.join(TROPoetree,tropoe_dir,'retrieval_output');
        f_pfx=tropoe_pid.lower()+'tropoe'+tropoe_dir+'.rt.2'
        dateindex=2
    return  sid_path, f_pfx, dateindex, suffix;

def get_gmlceil_info(project,TROPoetree,ceil_Sid,ceil_Pid):
    if 'wfip3' in project:
        suffix=".cdf"
        sid_path = os.path.join(TROPoetree,'CL_'+ceil_Pid.lower());
        f_pfx = 'ceil_'
        dateindex=1
    else:
        suffix=".cdf"
        sid_path = os.path.join(TROPoetree,'CL_'+ceil_Pid.lower());
        f_pfx = 'ceil_'
        dateindex=1
    return  sid_path, f_pfx, dateindex, suffix;

def read_var_fill_missing(dataset,varname,missing_var_to_mimic):
    if varname in dataset.variables:
        var = dataset.variables[varname]
        varvals = np.asarray(var);
        if 'missing_value' in var.ncattrs():
           varvals[varvals == var.missing_value] = np.nan
        elif '_FillValue' in var.ncattrs():
           varvals[varvals == var._FillValue] = np.nan
    else:
        varvals = np.full_like(missing_var_to_mimic, np.nan)

    return varvals

