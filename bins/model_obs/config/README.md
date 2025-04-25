## Configuration

When setting up a new project, there are several configuration files that need to be generated.  You can use the .default files as a guide to 
create project specific .py files

## config

#### config.py

__config.py__ specifies the disk location of the code, data, model data, and the TROPoe data.  This module
is imported by many of the programs to find the appropriate resources.  

CODEdir = ''   #  where the code resides, e.g., '/home/dmurray/scripts/WFIP3/modelobs'\
DATAdir = ''   #  where the data resides, e.g., '/Projects/WFIP3/modelobs'\
MODELdir = ''  # where the model data resides, e.g., '/psd3data/wfipdata/RAP_HRRR_columns'\
TROPoedir = '' # top level TROPoe data directory, e.g.,'/Projects/WFIP3/processed/realtime/TROPoe/'

Copy config.py.default to config.py and edit the entries as necessary for the project.

## config/parameter_files

#### forecast_param_list.py

__forecast_param_list.py__ specifie the number of timesteps available for each of the models used in the project. The model 
data in this case is the model profiles that are retrieved from GSL.

Copy forecast_param_list.py.default to forecast_param_list.py and edit the entries as necessary for the project.

#### model_param_list.py

__model_param_list.py__ specifies the models and stations used for processing the PSL generated ASCII model files. This is 
used by Model_Data/MDLdata_dwn.py to download the data and convert the ASCII to netCDF.  The generated netCDF files are
used for analyzing Disdrometer data since the GSL generated model profiles don't generally contain accumlated precipitation fields.

Copy model_param_list.py.default to model_param_list.py and edit the entries as necessary for the project.

#### model_station_list.py

__model_station_list.py__ defines the locations of the grid point locations for each station for the ASCII model files since they
don't contain that information.

Copy model_station_list.py.default to model_station_list.py and edit the entries as necessary for the project.

#### plotting_parameters.py

__plotting_parameters.py__ specifies many of the parameters that are used in the plotting of the images.  This includes max/min ranges for
various plots, averaging periods, plot locations.

Copy plotting_parameters.py.default to plotting_parameters.py and edit the entries as necessary for the project.

#### signal_param_list.py

__signal_param_list.py__ defines the stations of interest, the number of days to process, and which instruments (signals) are associated with
each station. 

Copy signal_param_list.py.default to signal_param_list.py and edit the entries as necessary for the project.

#### station_param_list.py

__station_param_list.py__ specifies information for each observing station, including location (lat/lon/alt), surface obs parameters, project id,
additional model grid location if using more than one for comparison, and the parameters to plot.

Copy station_param_list.py.default to station_param_list.py and edit the entries as necessary for the project.
