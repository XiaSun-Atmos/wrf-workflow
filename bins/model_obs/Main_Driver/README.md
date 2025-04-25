## Main_Driver

The Main_Driver directory contains the code for generating plots used in the model-obs comparison website for the project.  The code is broken out
into serveral subdirectories:

### driver_instr_functions

__driver_instr_functions__ contains subdirectories for each of the different instrument plotting code:

- __ESRL__  - code for reading the GSL generated model profile files
- __GMLceil__ - code for GML ceilometer data plots
- __PRFlidar__ - code for profiling lidar plots
- __PRFradar__ - code for PSL profiling radar plots
- __PSLdisd__ - code for PSL disdrometer plots
- __PSLmdl__ - code for reading the PSL generated model netCDF files
- __PSLmet__ - code for PSL surface met plots
- __TROPoe__ - code for TROPoe plots (both ASSIST and MWR)

### drivers2run

__drivers2run__ contains the code to generate plots for each instrument type.  

### general_functions

__general_functions__ conatins general code for analyzing and regridding data

### plotting_functions

__plotting_functions__ conatins code related to the plotting of data, including color maps, determining max/min for
plots, formatting/labeling.
