"""
Description
-----------

Utility methods for wrangling model data. This manipulating xarray DataArrays/Datasets,
manipulating latitude/longitude grids, special handling for model data
on isobaric grids (ie GFS), special handling for the RAP grid. 
This module is imported by, and used throughout, model_data.py.

"""

import logging
from pathlib import Path
import xarray as xr
import numpy as np
import column_extractor.calculations as calc
from column_extractor.file_io import file_exists

#### MANIPULATING/WRANGLING XARRAY DATAARRAYS/DATASETS

def rename_var(var_info: dict, var_ds: xr.Dataset) -> xr.Dataset:
    """Rename variable in xarray.Dataset

    We consider two renaming options:
        #. Specify a `rename_to` in variable_mapping.json (for custom renaming)
        #. If no `rename_to`, automatically rename variable to `shortName` provided
            in variable_mapping.json, to ensure consistency.
    Variable renaming is necessary to override cfgrib's automatic renaming of variables
    to CF conventions. Some of these variable names differ from GRIB naming conventions +
    what we prescribe in variable_mapping.json and can cause confusion.

    :param var_info: info for that variable, from variable_mapping.json
    :type var_info: dict
    :param var_ds: Dataset containing the variable we want to rename
    :type var_ds: xarray.Dataset
    :returns: Dataset containing the renamed variable
    :rtype: xarray.Dataset
    """

    # if variable_mapping.json contains instructions to do a custom variable renaming
    # (via rename_to), then rename the variable to that
    if 'rename_to' in var_info.keys():
        logging.info('    Renaming variable from %s to %s',
            list(var_ds.data_vars)[0], var_info['rename_to'])
        return var_ds.rename_vars({list(var_ds.data_vars)[0]: var_info['rename_to']})

    # if no explicit rename_to provided, rename the variable to the shortName prescribed
    # in variable_mapping.json. This clears up any confusion with cfgrib doing CF renaming
    if list(var_ds.data_vars)[0] != var_info['filter_keys']['shortName']:
        logging.info('    Renaming variable from %s to %s',
            list(var_ds.data_vars)[0], var_info['filter_keys']['shortName'])
        return var_ds.rename_vars(
            {list(var_ds.data_vars)[0]: var_info['filter_keys']['shortName']})

    return var_ds

def expand_vertical_dimension(coord_type: str, var_ds: xr.Dataset) -> xr.Dataset:
    """Expand the vertical coordinate dimension, if needed

    Variables with a single vertical level (ie surface, cloudBase, nominalTop, etc.) are
    forced into 2D lat/lon arrays. These can't be easily merged into a single xarray.Dataset
    with like vertical coordinates (since it's like the vertical coordinate doesn't exist).
    So, we expand along this dimension to convert these 2D lat/lon arrays into 3D
    vertical/lat/lon arrays of size [1, lat, lon]. These now 3D variables can be merged into
    a single Xarray.Dataset with like vertical coordinates.

    :param coord_type: vertical coordinate type (ie hybrid, surface, etc.)
    :type coord_type: str
    :param var_ds: Dataset containing the variable we may / may not want to expand dimensions
    :type var_ds: xarray.Dataset
    :returns: Dataset containing the variable with (potentially) expanded vertical coordinate
    :rtype: xarray.Dataset
    """

    # for vertical coordinates with only one level (ie surface, low cloud layer, etc)
    # expand the vertical dimension so that they can be effectively merged
    # along the vertical dimension (ie from 2D to 3D)
    if var_ds.coords[coord_type].size == 1:
        logging.info('    Vertical coordinate is size 1, so expanding to from 2D to 3D')
        return var_ds.expand_dims(coord_type)

    # if vertical coordinate contains multiple levels (ie hybrid/isobaric/height levels),
    # then no expansion needed
    return var_ds

def keep_specific_vertical_levels(
    coord_type: str,
    var_info: dict,
    var_ds: xr.Dataset) -> xr.Dataset:
    """Only retain specific vertical levels, based on `keep_levels` or `drop_levels`

    Option to specify levels to keep OR levels to drop. This is done using xr.Dataset.sel()
    and xr.Dataset.drop_sel(), respectively. Note that values in `keep_levels`/`drop_levels`
    in variable_mapping.json must exactly match values associated with the current vertical
    coordinate + must match its data type and units. EG, one would pass
    `keep_levels = [1000.,975.,950.,925.,900.]` to retain the pressure levels between 1000 hPa
    and 900 hPa and drop all others. One would pass `drop_levels = [1000.,975.,950.,925.]`
    to drop all levels below 900 hPa.

    :param coord_type: vertical coordinate type (ie hybrid, surface, etc.)
    :type coord_type: str
    :param var_info: info for that variable, from variable_mapping.json
    :type var_info: dict
    :param var_ds: Dataset containing the variable we want to rename
    :type var_ds: xarray.Dataset
    :returns: Dataset containing the renamed variable
    :rtype: xarray.Dataset
    """

    # if `keep_levels` is defined in variable_mapping.json for given variable,
    # then select only those levels
    if 'keep_levels' in var_info.keys():
        logging.info('    Only retaining the following vertical levels: %s',
            (', ').join([str(level) for level in var_info['keep_levels']]))
        return var_ds.sel({coord_type:var_info['keep_levels']})

    # if `drop_levels` is defined in variable_mapping.json for given variable,
    # then drop those levels and keep all others
    if 'drop_levels' in var_info.keys():
        logging.info('    Dropping the following %s levels: %s',
            coord_type,
            (', ').join([str(level) for level in var_info['drop_levels']]))
        return var_ds.drop_sel({coord_type:var_info['drop_levels']})

    return var_ds

def merge_vars_across_files(vars_across_files: list[xr.Dataset]) -> xr.Dataset:
    """Merge variable Datasets from different files into single variable Dataset

    Some variables are duplicated in different file types. Some variables have certain
    levels in one file and other levels in a different file. Here we use Xarray's
    merge function to synthesize these variables across files by reconciling duplicates,
    combining various vertical levels, etc. This approach takes the union of the datasets.

    :param vars_across_files: list of variable Datasets from different files
    :type vars_across_files: list
    :returns: merged Dataset containing single variable reconciled across multiple files
    :rtype: xarray.Dataset
    """

    logging.info('    Merging the above into a single Xarray Dataset')

    # take the list containing all variables for the given vertical coordinate
    # and use Xarray's combine_by_coords to combine these variables into a single Xarray Dataset
    return xr.merge(vars_across_files, combine_attrs='drop_conflicts')

def combine_vars_with_like_coords(vars_for_coord: list[xr.Dataset]) -> xr.Dataset:
    """Combine individual variable Datasets with like vertical coordinates into single Dataset

    As mentioned above, it's necessary to read variables one at a time, to avoid conflicts between
    similar variable names with different vertical coordinates, step types, etc. We
    want to combine these individual variable Datasets with like vertical coordinates into a
    single Dataset. This method uses Xarray's combine_by_coords function to accomplish this.

    :param vars_for_coord: list of variable Datasets with like vertical coordinates
    :type vars_for_coord: list
    :returns: combined Dataset containing all variables with like vertical coordinates from list
    :rtype: xarray.Dataset
    """

    logging.info('Combining the above into a single Xarray Dataset\n')

    # take the list containing all variables for the given vertical coordinate
    # and use Xarray's combine_by_coords to combine these variables into a single Xarray Dataset
    return xr.combine_by_coords(vars_for_coord, combine_attrs='drop_conflicts') # type: ignore

def _force_latitude_longitude_2d(merged_vars_by_coord_ds: xr.Dataset) -> xr.Dataset:
    """Force regular latitude/longitude grid to have 2D latitude/longitude arrays

    :param merged_vars_by_coord_ds: Dataset containing all variables for vertical coordinate
    :type merged_vars_by_coord_ds: xarray.Dataset
    :returns: Dataset with 2D latitude/longitudes (might not be necessary)
    :rtype: xarray.Dataset
    """

    # create 2D arrays of latitude, longitude
    longitude, latitude = np.meshgrid(
        merged_vars_by_coord_ds['longitude'].data, merged_vars_by_coord_ds['latitude'].data)

    # replace current 1D latitude/longitude arrays with the 2D versions
    for dim, values in {'latitude':latitude,'longitude':longitude}.items():
        merged_vars_by_coord_ds[dim] = xr.DataArray(
            values,
            dims=['y','x'],
            attrs=merged_vars_by_coord_ds[dim].attrs)

    # swap latitude/longitude dimensions for y,x
    return merged_vars_by_coord_ds.swap_dims({'latitude':'y','longitude':'x'})

def construct_step_range(filter_keys: dict[str, str], forecast_step: int) -> dict[str, str]:
    """Construct dynamic step range filter key to be used when reading in accumulated variables

    Some variables have an accumulated step type (`accum`), meaning they're accumulated
    over specific time ranges (instead of being instantaneous). EG: accumulated precipitation.
    Accumulated precipitation can come in two forms, depending on the forecast hour:
        - Precipitation accumulated since initialization (stepRange `0-N`, EG: `0-3`)
        - Precipitation accumulated over the last forecast hour (stepRange `N-1-N`, EG: `2-3`)
    where N is the current forecast hour. Here, we dynamically create the step range to filter
    against, based on the current forecast step and what is prescribed in the variable's
    filter_keys section in variable_mapping.json.

    :param filter_keys: original filter keys, containing placeholder step range filter key
    :type filter_keys: dict
    :param forecast_step: current forecast hour step (ie 0, 1, 2, 3, etc.)
    :type forecast_step: int
    :returns: filter keys with amended step range filter key (if needed)
    :rtype: dict
    """

    # simple case: forecast step is 0, meaning stepRange isn't a valid key
    # (but accumulated precipitation variable still exists). Drop from filter_keys:
    if forecast_step == 0:
        # drop stepRange from filter keys
        del filter_keys['stepRange']
    # for any other forecast hour, construct step range based on placeholder format
    elif forecast_step > 0:
        # could be 0-N: variable accumulated from initialization to current forecast step
        if filter_keys['stepRange'] == '0-N':
            # replace N with current forecast step
            filter_keys['stepRange'] = filter_keys['stepRange'].replace('N',str(forecast_step))
        # could be (N-1)-N: variable accumulated from previous forecast step to current
        elif filter_keys['stepRange'] == '(N-1)-N':
            # do some math to create (N-1) - N in string format
            filter_keys['stepRange'] = f'{forecast_step-1}-{forecast_step}'

    return filter_keys

#### LATITUDE/LONGITUDE GRID WRANGLING

def prepare_latitude_longitude(merged_vars_by_coord_ds: xr.Dataset) -> xr.Dataset:
    """Prepare latitude/longitudes for analysis

    This includes:
    - For regular latitude/longitude grids that are respectively 1D, create 2D
        latitude/longitude arrays and add to Dataset, changing dimensions from
        latitude/longitude to y/x (necessary for GFS)
    - Recast latitude/longitudes from 0-180 / 0-360 to -90 to 90 / -180 to 180, respectively)
        (necessary for HRRR/GFS longitudes)

    :param merged_vars_by_coord_ds: Dataset containing all variables for vertical coordinate
    :type merged_vars_by_coord_ds: xarray.Dataset
    :returns: Dataset with prepped latitude/longitudes (might not be necessary)
    :rtype: xarray.Dataset
    """

    # if Dataset dimensions include latitude/longitude (as opposed to y/x), then grid type
    # must be regular latitude/longitude, which means we need to force the latitude/longitude
    # arrays to be 2D instead of 1D
    if any(dim in merged_vars_by_coord_ds.dims for dim in ['latitude','longitude']):
        merged_vars_by_coord_ds = _force_latitude_longitude_2d(merged_vars_by_coord_ds)

    # recast longitudes from 0 to 360 to -180 to 180, if needed
    merged_vars_by_coord_ds['longitude'] = merged_vars_by_coord_ds['longitude'].where(
        merged_vars_by_coord_ds['longitude'] <= 180.,
        merged_vars_by_coord_ds['longitude'] - 360.)

    # recast latitudes from 0 to 180 to -90 to 90, if needed
    merged_vars_by_coord_ds['latitude'] = merged_vars_by_coord_ds['latitude'].where(
        merged_vars_by_coord_ds['latitude'] <= 90.,
        merged_vars_by_coord_ds['latitude'] - 180.)

    return merged_vars_by_coord_ds

def add_rot_grid(merged_vars_by_coord_ds: xr.Dataset, model: str) -> xr.Dataset:
    """Calculate and add rotation grid to coordinates; behavior depends on model

    We need to convert model winds from grid-relative to earth-relative, which first requires
    computing the grid rotation (in radians) at each latitude/longitude coordinate, if any.
    This varies depending on model grid projection:
    - HRRR: compute for Lambert Conformal via https://rapidrefresh.noaa.gov/faq/HRRR.faq.html
    - GFS: compute for regular latitude/longitude grid, which requires no transformation, so
        return array of 0s which acts as the zero matrix.
    - RAP: we already have the rotation grid and apply this in a different method. So return
        original data

    :param merged_vars_by_coord_ds: Dataset containing all variables for vertical coordinate
    :type merged_vars_by_coord_ds: xarray.Dataset
    :returns: Dataset with model rotation grid added to coordinates
    :rtype: xarray.Dataset
    """

    # calculate HRRR rotation grid
    if model == 'hrrr':
        rotation = calc.get_hrrr_rotation_grid(merged_vars_by_coord_ds['longitude'])
    # GFS rotation grid is zero matrix, since regular latitude-longitude grid
    elif model == 'gfs':
        rotation = xr.zeros_like(merged_vars_by_coord_ds['longitude']).rename('rotation')
        rotation.assign_attrs({
            'standard_name':'rotation',
            'long_name':'GFS rotation grid (zero matrix; since regular lat-lon)'})
    # RAP rotation grid already exists from previous method, pass
    elif model == 'rap':
        return merged_vars_by_coord_ds
    elif model == 'wfip3km' or model == 'wfip1km':
        rotation = calc.get_wfip3_rotation_grid(merged_vars_by_coord_ds['longitude'])

    return merged_vars_by_coord_ds.assign_coords({'rotation':rotation})

#### SPECIAL NEEDS FOR CERTAIN MODELS (IE RAP GRID) AND VERTICAL COORDINATES (IE PRESSURE LEVELS)

def prep_isobaric_data(isobaric_dataset: xr.Dataset) -> xr.Dataset:
    """Prep data on isobaric levels (create pressure variable, flip order)

    :param isobaric_data: isobaric data to be prepped
    :type isobaric_data: xarray.Dataset
    :returns: isobaric data with pressure variable and reversed vertical coordinate
    :rtype: xarray.Dataset
    """

    logging.info('Prepping isobaric data (create pressure variable, flip order)...')

    logging.info('--- Creating pressure variable (same profile at all gridpoints)')

    # create a pressure variable (which is technically the same vertical profile
    # at each latitude/longitude point - but is important when interpolating to heights)
    isobaric_dataset['pres'] = isobaric_dataset['isobaricInhPa'].expand_dims(
        {'y':isobaric_dataset.sizes['y'],'x':isobaric_dataset.sizes['x']})

    # reorder dimensions of pressure variable to isobaricInhPa, y, x
    isobaric_dataset['pres'] = isobaric_dataset['pres'].transpose('isobaricInhPa','y','x')

    logging.info('--- Reversing pressure level order to highest -> lowest\n')

    # reindex the Dataset such that pressure levels are from highest (ie lowest altitude)
    # to lowest (ie highest altitude) so that we don't have to reverse them later
    isobaric_dataset = isobaric_dataset.reindex(
        {'isobaricInhPa':isobaric_dataset['isobaricInhPa'][::-1]})

    return isobaric_dataset

def get_rap_lat_lon_grid(rap_lat_lon_grid_path: Path) -> xr.Dataset:
    """If RAP grid netCDF file exists, read it in

    cfgrib/eccodes/xarray can't handle RAP's outer domain's Arakawa-E grid and returns a 1D
    lat/lon array. We need a standalone grid to supplement, which is where this comes in.
    File from https://noaa-ufs-srw-pds.s3.amazonaws.com/index.html#fix/fix_am/

    :param rap_lat_lon_grid_path: path to RAP lat/lon grid file (should be in conf directory)
    :type rap_lat_lon_grid_path: pathlib.Path
    :returns: Dataset cotaining RAP lat/lon/rot grids
    :rtype: xarray.Dataset
    """

    logging.info('Reading %s...\n', rap_lat_lon_grid_path.name)

    # check if file exists
    file_exists(rap_lat_lon_grid_path)

    # read netcdf file containing RAP lat lon grid
    rap_lat_lon_grid = xr.open_dataset(rap_lat_lon_grid_path, engine='netcdf4')

    # rename lat/lon grid dimensions from ny,nx to y,x
    #to work with default xarray conventions
    rap_lat_lon_grid = rap_lat_lon_grid.rename_dims({'ny':'y','nx':'x'})

    return rap_lat_lon_grid

def add_rap_lat_lon_grid(
    rap_lat_lon_grid: xr.Dataset,
    coord_type: str,
    var_: xr.DataArray) -> xr.DataArray:
    """Add RAP lat/lon grid to variable

    As previously mentioned, the RAP outer domain's grid isn't recognized by cfgrib/eccodes.
    As such, we need to add it after the fact. Here, we reshape the variable from
    [vertical coordinate, values] to [vertical coordinate, latitude, longitude] using the
    dimensions of the lat/lon grid that we read in. Then, we assign new latitude/longitude
    coordinates using the lat/lon grid and return a DataArray with these new coordinates.

    :param rap_lat_lon_grid: Dataset containing the RAP lat/lon grid
    :type rap_lat_lon_grid: xarray.Dataset
    :param coord_type: vertical coordinate type (ie hybrid, surface, etc.)
    :type coord_type: str
    :param var_: the variable we want to add the RAP lat/lon grid to
    :type var_: xarray.Dataset
    :returns: DataArray containing variable with newly added lat/lon coordinates
    :rtype: xarray.DataArray
    """

    # reshape variable from 2D [vertical coord, values] to 3D [vertical coord, y, x]
    # using the shape of the latitude/longitude grids
    reshaped_var = var_.data.reshape(
        tuple([var_[coord_type].size]) + rap_lat_lon_grid['gridlat'].shape)

    # define new coordinates of:
    # vertical coord, longitude, latitude, and rotation grids
    new_coords = {
        coord_type:var_[coord_type],
        'longitude':rap_lat_lon_grid['gridlon'],
        'latitude':rap_lat_lon_grid['gridlat'],
        'rotation':rap_lat_lon_grid['gridrot']}

    # now convert this reshaped numpy array to an Xarray dataArray
    # assign dimensions as [vertical coord, y, x] and coordinates as
    # vertical coordinate, longitude (longitude grid), latitude (latitude grid), rotation
    return xr.DataArray(reshaped_var, dims=[coord_type,'y','x'], coords=new_coords)
