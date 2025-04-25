

# Overview

This application currently runs in real-time in GSL's production AWS environment. For details on this
infrastructure, see
[dsg-cloud-pipelines repository](https://github.com/NOAA-GSL/dsg-cloud-pipelines/tree/main/column_extraction_workflow).

The Python code, configuration files, and supplementary data needed to run the application
are bundled into a Docker container image.

The "cookbook" for building this container
image is outlined in the [Dockerfile](../Dockerfile).

Container image versions are directly tied to tags in the repository. The 5 most recent container
images are stored in Harbor (dev) and AWS ECR (dev2 and prod). They can also be built locally. For
more info, see [build instructions](build_and_deploy.md).

For a history of versions, see
[tags](https://github.com/NOAA-GSL/column-extractor/tags) and
[releases](https://github.com/NOAA-GSL/column-extractor/releases).

## Repository layout

Source code and tests sit under [src](../src).

Configuration files and supplementary data sit under [conf](../conf).

The [Dockerfile](../Dockerfile) outlines how the container image is built and what's contained therein.

[requirements.yml](../requirements.yml) outlines the Python dependencies required.

Miscellaneous scripts for building/pushing new versions of the container image, running our tests,
etc. sit under [integration](../integration).

Documentation sits under [doc](../doc).

## Application workflow

The application performs the following tasks:

### Pre-process the data

1. For each vertical coordinate type defined in `variable_mapping.json`:
   extract each variable of interest, rename it (if necessary), reshape
   it (if necessary), filter vertical levels (if necessary), and append
   its Xarray Dataset to a list.
2. Merge all variables for a given vertical coordinate type into a single
   Xarray Dataset; add to dictionary with vertical coordinate type as key.
3. Add latitude/longitude grid from file to each Xarray Dataset if isn't
   recognized by cfgrib/eccodes (ie RAP)
4. Recast latitude/longitude grids as needed (ie 0-360deg to -180-180deg)
5. Compute and/or add rotation grid (relevant for converting grid-relative
   winds to earth-relative)
6. If pressure level data, do additional pre-processing

### Extract data at closest model gridpoint to each site of interest

1. Identify the closest model gridpoint to each site of interest using
   the minimum great circle distance calculated via the Haversine function.
2. Subset the data to retain only an NxN neighborhood surrounding the
   closest model gridpoint to each site of interest. NxN is determined
   by `neighbors` in `column_extractor.cfg` (ie `neighbors = 1` yields
   a 3x3 neighborhood, `neighbors = 2` yields a 5x5 neighborhood, etc.)

### Post-process the subsetted data

1. (if hybrid levels) derive elevation above MSL, which Dave defines as
   the geopotential height at the first hybrid level
2. Scale geopotential height such that the first level is 5m and all other
   levels are relative to that (another Dave convention)
3. Derive mixing ratio at 2m and on all hybrid levels from specific
   humidity
4. Derive virtual temperature on all hybrid levels from temperature and
   mixing ratio
5. Derive geometric vertical velocity (ie m/s) on all hybrid levels
   from pressure vertical velocity (ie Pa/s), pressure, and temperature
6. Convert grid-relative winds at 2/10/80m and on all hybrid levels to
   earth-relative winds
7. Interpolate 3D variables onto standardized height grid (defined by
   `new_heights_m` in `column_extractor.cfg`)
8. Perform any unit conversions, as defined by `convert_to` in
   the `post-process` sections of `variable_mapping.json`

### Interpolate the subsetted/post-processed data to each site of interest

1. Loop through each variable and interpolate each NxN neighborhood to the
   corresponding site of interest using Scipy's `griddata()`

### Write data to netCDF file

1. Generate temporary CDL file from CDL template (fill in placeholders)
2. Generate output netCDF file from this temporary CDL file
3. Populate the output netCDF file with the closest model gridpoint data
   and interpolated data for each variable. `netcdf_mapping` in the
   `post-process` sections of `variable_mapping.json` dictates how each
   variable is mapped to a variable in the output netCDF file

## How do I build the container image? How do I deploy a new version?

See [build instructions](build_and_deploy.md).

## How do I run the application?

See [run instructions](run.md).

## How do I add support for a new model, variable, or site of interest?

See [configuration instructions](configure/README.md).

## How do I debug when things go wrong?

See [debug instructions](debug.md).
