# Overview

All configuration files associated with the application sit under
[conf](../../conf).

The configuration changes you intend to make dictate which files
need modified/added and to what extent.

There are three common configuration changes:

- To add a site for one or more models, see [doc](add_a_site.md)
- To add a variable for one or more models, see [doc](add_a_variable.md)
- To add a new model, see [doc](add_a_model.md)

Each configuration file is discussed below.

## [`column_extractor.cfg`](../../conf/column_extractor.cfg)

Contains a section for each model that defines:

- `vertical_coord`: the vertical coordinate type (ie `hybrid`, `isobaricInhPa`, etc.)
- `new_heights_m`: the height levels to interpolate the 3D data to
  (10m, 20m, 40m, ..., 15km)
- `neighbors`: how many neighbors to use when performing billinear interpolation (typically 1)
- `lat_lon_grid`: (optional) the path to the latitude/longitude grid for this model,
  if it cannot be derived from the input files themselves

The RAP section looks like:

```python
[rap]
vertical_coord = hybrid
new_heights_m = 10,20,40,60,80,100,120,140,160,180,200,225,250,275,300,350,400,450,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4200,4400,4600,4800,5000,5500,6000,6500,7000,8000,9000,10000,11000,12000,13000,14000,15000
lat_lon_grid = ./conf/rap_arakawa_grid_3.32769.nc
neighbors = 1
```

> [!NOTE]
> RAP is the only model we have a latitude/longitude grid file for,
> because its full North American domain isn't recognized by
> cfgrib/eccodes and cannot be computed on the fly.

## [`profiling_sites.json`](../../conf/profiling_sites.json)

Contains the names, latitudes, and longitudes for each site 
we wish to extract model data at. There's a section for each model.

The file has structure:

```json
{
    "model": {
        "site name": [
            "site latitude (degrees N)",
            "site longitude (degrees E)"
        ]
    }
}
```

For example:

```json
{
    "hrrr": {
        "ARM sgpC1 site (Lamont, OK)": [
            36.6053,
            -97.4857
        ]
    }
}
```

## [`variable_mapping.json`](../../conf/variable_mapping.json)

This file contains `pre-process` and `post-process` sections for
each model.

The `pre-process` section defines parameters used during the pre-processing step:

- Every variable to extract from input file(s), defined by its
  `shortName`, `stepType`, and the vertical coordinate type it sits
  under (ie `hybrid`, `surface`, etc.)
- Which input file(s) to extract it from (via `file_id`)
- (optional) what to rename the variable to (via `rename_to`)
- (optional) which levels to drop (via `drop_levels`) or keep (via `keep_levels`)

The `post-process` section defines parameters used during the post-processing and output writing
steps:

- How a variable maps to its counterpart(s) in the output netCDF file (via `netcdf_mapping`)
- Any unit conversions to perform before writing to output netCDF file (via `convert_to`)

### Extracting specific variables from GRIB files

Since we don't need all variables in a GRIB file, we extract each variable of interest
individually. This is the main purpose of the `pre-process` section.

To ensure that only one variable is targeted at a time, we leverage
Xarray/cfgrib's `filter_keys` parameter. The main `filter_keys` used are:

- `shortName`: the "name" of the variable
- `stepType`: the type of time step (ie `instant`, `accum`, `avg`, etc.). This will almost always
  be `instant`, but can be `accum` for accumulated precipitation variables and `avg` for variables
  that are averaged over the previous forecast hours.
- `typeOfLevel`: the vertical coordinate type (ie `hybrid`, `surface`, `cloudBase`, etc.)
- (occasionally) `stepRange`: used with `stepType` `accum` to define which forecast hour range
  we want. GRIB files often have accumulated precipitation from FH00 to FHNN and FHN-1 to FHNN. We
  typically only want FHN-1 to FHNN. We use `"stepRange": "(N-1)-N"` and let the application
  dynamically populate the correct step range.

> [!NOTE]
> The taxonomy for all of the above is defined by cfgrib/eccodes and follows CF Conventions.
> They have no correlation to those listed in NCEP `.idx` files. Some are easily
> "translatable", others aren't.

## [`columns_output_{hrrr,rap,gfs}_template.cdl`](../../conf/columns_output_hrrr_template.cdl)

These CDL files are templates used to generate output netCDF files for each forecast hour.
There's one for each model, since certain variables only exist for certain models.

They define:

- The dimensionality of the output file
- Each variable and its metadata (long names, units, fill values)
- Global metadata (source file, model, code version, etc.)

For example:

```
float lwc0(forecast_hour, sites, height_bin) ;
        lwc0:long_name = "Liquid water mixing ratio at the closest model grid location" ;
        lwc0:units = "kg/kg" ;
        lwc0:_FillValue = 9.96921e+36f ;
float lwc1(forecast_hour, sites, height_bin) ;
        lwc1:long_name = "Liquid water mixing ratio calculated by bilinear interpolation via Scipy's griddata()" ;
        lwc1:units = "kg/kg" ;
        lwc1:_FillValue = 9.96921e+36f ;
```

Defines liquid water mixing ratio at the closest model grid location (`lwc0`) and interpolated to
each site (`lwc1`). The suffix `0` denotes closest model grid location and `1` denotes interpolated
to each site.

The applicable `post-process` section of `variable_mapping.json` will map the variable's `shortName` to its name in the output netCDF file:

```json
{
    "var_name": "clwmr",
    "netcdf_mapping": ["lwc"]
}
```

> [!NOTE]
> The variable names and `0`/`1` suffixes are legacy conventions from Dave's output. The application
> appends the suffixes automatically when writing to netCDF.

In the dimensions and global attributes, you'll see several instances of all caps. EG:

```
dimensions:
        forecast_hour = unlimited ;
        sites = N_SITES ;
        height_bin = N_HEIGHT_BINS ;
        depth_bin = N_DEPTH_BINS ;
```
and
```
// global attributes:
                :title = "Vertical profiles / point values at profiling sites of interest" ;
                :institution = "NOAA Global Systems Laboratory (GSL)" ;
                :contact = "dsg.its.gsl@noaa.gov, dave.turner@noaa.gov" ;
                :source = SOURCE_FILE ;
                :model = "HRRR V4" ;
	            :model_initialization = MODEL_INIT ;
                :version = CODE_VERSION ;
                :references = "https://github.com/NOAA-GSL/column-extractor" ;
                :Conventions = "CF-1.10" ;
                :date_file_created = FILE_CREATION_TIME ;
```

These placeholders are dynamically populated when a temporary CDL file is
created from the template CDL. The output netCDF file is then generated
from this temporary CDL file. These serve two purposes:

- We don't need to change the CDL files if we add/remove a profiling site or change the
  vertical height grid we interpolate to (thus changing dimension)
- Each output file has unique creation time, input file(s), model initialization, and
  code version metadata

## Custom GRIB definitions under [`conf/custom_definitions`](../../conf/custom_definitions)

Depending on the model output, eccodes version, and GRIB tables version some variables may
not be decoded. Some or all of the attributes of this variable will be listed as `unknown`.
Output from `grib_ls` might look like:

```
edition    centre    date    dataType    gridType    stepRange    typeOfLevel    level    shortName    packingType  
2    kwbc    20240603    fc    lambert    1    surface    0    unknown    grid_complex_spatial_differencing 
```

Note the `unknown` for `shortName`.

In order to correct this, you'll need to add the correct GRIB codes under
`conf/custom_definitions/grib2/localConcepts/kwbc` (`kwbc` may change depending on `centre`). For directions
on this, see [doc](../debug.md#fixing-unknown-grib-variables).
