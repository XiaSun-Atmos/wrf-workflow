# PSL obs decoders and Utilities

This is a collection of decoders and utilities for reading in the PSL realtime observation (/psd2data/obs/realtime) ASCII files and outputting netCDF files.  

The following decoders are available:

- psldisd_to_netcdf.py - convert Disdrometer data (/psd2data/obs/realtime/DisdrometerParsivel/Stats)
- pslmdlmet_to_netcdf.py - convert PSL extracted model data 
- pslrwp_to_netcdf.py    - convert Radar Wind Profiler data (/psd2data/obs/realtime/Radar[449|915]/WwWind\*)
- pslsfcmet_to_netcdf.py - convert surface datalogger data (/psd2data/obs/realtime/CsiDatalogger/SurfaceMet)

Variable config files:
- pslobs_vars_dict.py - dictionary of obs variable parameters
- mdlobs_vars_dict.py - dictionary of model variable parameters

Other utilities:
- pslutils.py
- psldateutils.py

