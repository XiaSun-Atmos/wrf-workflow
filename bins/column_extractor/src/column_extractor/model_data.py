"""
Description
-----------

Class ModelData handles all aspects of model data handling: extracting all variables from file
(as prescribed by variable_mapping.json), merging variables with like vertical coordinates into
single Xarray Dataset (doing all necessary coordinate massaging), any model-specific preprocessing
(adding in RAP grid, etc.), adding each vertical coordinate Xarray Dataset into coordinate-wise
dictionary, as well as all post-processing (derivations, interpolation, etc.).

"""

import logging
from pathlib import Path
import xarray as xr
import numpy as np
from column_extractor.proc_config import ProcConfig
import column_extractor.model_utils as utils
import column_extractor.calculations as calc
from column_extractor import file_io

class ModelData(ProcConfig):
    """Handles all model data processing (preprocessing, variable calculations, column extraction)

    Specific processing nuances depend on the model the input file is associated with (IE need to
    add grid for RAP, need to fix latitude/longitude for HRRR, some model fields vary, etc.)

    Attributes:
        proc_config (class): ProcConfig class instance; contains all processing configuration info

    Methods:
        preprocess(): Handles initial raw data processing/massaging
        extract_closest_gridpoint_data(): Extracts data at closest model gridpoint to each site
        post_process(): Derives additional fields, performs unit conversions, etc.
        interpolate_to_profiling_sites(): Interpolates model data to each profiling site
        write_output_to_netcdf(): Writes output to netCDF file, following CDL conventions
    Usage:
        # Create an instance of ModelData
        model_data = ModelData(proc_config)

        # read in model file and pre-process the data
        model_data.preprocess()

        # extract model data at closest model gridpoint to each profiling site
        model_data.extract_closest_gridpoint_data(preprocessed_data)

        # post-process the data
        model_data.post_process(closest_gridpoint_data)

        # interpolate post-processed data to each profiling site
        model_data.interpolate_to_profiling_sites(post_processed_data)

        # write output data to netCDF file, following CDL conventions
        model_data.write_output_to_netcdf(post_processed_data, interpolated_data)

        # Access the attributes and methods of the instance
        print(model_data.<attribute>)
    """

    def __init__(self) -> None:
        """Initialize the ModelData instance and add proc_config information as attribute

        :param proc_config: path to RAP lat/lon grid file (should be in conf directory)
        :type proc_config: dict
        :returns: {None}
        :rtype: {None}
        """

        # inherit model-specific configuration parameters for use throughout ModelData
        ProcConfig.__init__(self)

    def _drop_var_from_post_process(self, coord_type: str, var_info: dict) -> None:
        """Drop variable from post-process section of variable mapping, so we don't reference it

        Some variables aren't available at FH00 vs. other forecast hours (and vice versa). If we
        attempt to extract a variable that isn't present, remove it from the post-process
        section of variable mapping, so that we don't run into issues in the rest of the code.
        That way, nothing will be written to that variable in the output file and it'll be filled
        with its _FillValue

        :param coord_type: vertical coordinate type (ie hybrid, surface, etc.)
        :type coord_type: str
        :param var_info: variable to drop's metadata: filter keys, rename_to, etc.
        :type var_info: dict
        :returns: {None}
        :rtype: {None}
        """

        # determine which variable to drop
        # could be new name, as defined by rename_to
        if 'rename_to' in var_info:
            var_to_drop =  var_info['rename_to']
        # if no rename_to specified, then we use standard shortName
        else:
            var_to_drop = var_info['filter_keys']['shortName']

        # log which variable we're going to remove
        logging.warning("%s couldn't be found in file! Removing %s from post-process...",
            var_info['filter_keys']['shortName'],
            var_to_drop)

        # loop through each variable in post-process section for given coordinate type
        for var_ in self.model_var_map['post-process'][coord_type]:
            # if variable's name matches the variable to drop, drop it
            if var_['var_name'] == var_to_drop:
                self.model_var_map['post-process'][coord_type].remove(var_)

    def _extract_var_from_files(self, filter_keys: dict[str, str], file_id: list) -> xr.Dataset:
        """Extract variable from each input file(s), if applicable

        We may receive one concatenated GRIB file as input or several GRIB files that,
        together, contain all the fields we're interested in. Here, we loop through each
        file and extract the variable in question from the file it resides in. If the variable
        exists across multiple files, then we merge these individual Datasets.

        :param filter_keys: set of keys representing a single, unique variable in a file
        :type filter_keys: dict
        :param file_id: list containing file identifiers corresponding to file variable is in
        :type file_id: list
        :returns: Dataset containing the variables merged from multiple files
        :rtype: xarray.Dataset
        """

        # declare empty list to add variable Dataset from each input file
        var_ds_across_files = []

        # loop through each input file
        for file in self.input_output['local_paths']:
            # if variable is in this input file (or no target file specified), retrieve
            if not file_id or any(file_type in str(file) for file_type in file_id):
                logging.info('From file: %s', file)
                # extract variable from current file
                var_ds_per_file = file_io.read_grib(file, filter_keys)
                # add to list for merging
                var_ds_across_files.append(var_ds_per_file)

        # if merge list of same variable from different files into single Dataset
        if len(var_ds_across_files) > 1:
            return utils.merge_vars_across_files(var_ds_across_files)

        # otherwise, just return first (and only) entry in list
        return var_ds_across_files[0]

    def _extract_vars_for_vertical_coord(
        self,
        coord_type: str,
        vars_to_extract: list[dict]) -> tuple[dict, xr.Dataset]:
        """For each variable in vertical coordinate type: read, rename, expand, etc.

        We want a range of variables for each vertical coordinate type. Here, we extract each
        variable from the file(s), rename it if needed, expand its vertical coordinate dimension
        if needed, filter out unwanted vertical levels if needed, and add it to a dictionary with
        key being variable name. Method then returns this dictionary and the merged Dataset
        containing all of the variables.

        :param coord_type: vertical coordinate type (ie hybrid, surface, etc.)
        :type coord_type: str
        :param vars_to_extract: the variables we're interested in for given `coord_type`
        :type vars_to_extract: list
        :returns: tuple containing dict with each variable and merged Dataset of all variables
        :rtype: (dict, xarray.Dataset)
        """

        logging.info('Extracting variables for vertical coordinate type: %s', coord_type)

        # declare empty list for the given vertical coordinate's variables
        vars_for_coord = []

        # declare empty dict to hold each variable's attributes
        var_attrs = {}

        # declare forecast step as None, to be filled by first variable extraction
        forecast_step = None

        # loop through each variable for this vertical coordinate; extract from file
        for var_info in vars_to_extract:
            # create filter_keys dictionary to pass to _read_file()
            filter_keys = {'typeOfLevel':coord_type} | var_info['filter_keys']
            # if stepRange placeholder is specified in filter_keys, dynamically fill it
            # (but make sure we've already set forecast_step, since we need it)
            if 'stepRange' in filter_keys and forecast_step is not None:
                filter_keys = utils.construct_step_range(filter_keys, forecast_step)
            # extract variable from file(s), merge into single variable Dataset
            var_ds = self._extract_var_from_files(filter_keys, var_info['file_id'])
            # if we don't have the forecast hour yet, get it and set it
            # (for use in any variables with accumulated stepTypes and specific stepRanges)
            if forecast_step is None and 'step' in var_ds.variables:
                forecast_step = calc.forecast_hour_from_timestep(var_ds['step'])
            # verify that variable was retrieved from file successfully
            if len(var_ds.data_vars) == 0:
                # if not, drop it from post-process section (so we don't reference it)
                self._drop_var_from_post_process(coord_type, var_info)
            else:
                # handle any variable renaming that may be needed
                var_ds = utils.rename_var(var_info, var_ds)
                # for vertically 1D fields (ie surface), expand the vertical dimension as 1
                var_ds = utils.expand_vertical_dimension(coord_type, var_ds)
                # filter out specific vertical levels, if needed
                var_ds = utils.keep_specific_vertical_levels(coord_type, var_info, var_ds)
                # add each variable to list (for merging)
                vars_for_coord.append(var_ds)
                # add each variable's attributes to dictionary (for units/names/etc.)
                var_attrs[list(var_ds.data_vars)[0]] = (var_ds[list(var_ds.data_vars)[0]].attrs |
                    {'vertical_coordinate':coord_type})

        # combine variables with like vertical coordinate into single Xarray Dataset
        return var_attrs, utils.combine_vars_with_like_coords(vars_for_coord)

    def preprocess(self) -> dict[str, xr.Dataset]:
        """Pre-process the raw model data based on the provided ProcConfig instance

        This includes:
        - Read the file in
        - Extract variables for each vertical coordinate
        - Merge variables into single Xarray Dataset for each vertical coordinate
        - Adds vertical coordinate Xarray Dataset to vertical coord-wise dictionary.

        :param {None}: {None}
        :returns: Vertical coord-wise dictionary containing Xarray Datasets of preprocessed data
        :rtype: dict
        """

        # raise cfgrib logger threshold from WARNING to ERROR
        # this suppresses warnings associated with eccodes not handling the RAP grid
        logging.getLogger('cfgrib').setLevel(logging.ERROR)

        logging.info(
            'Preprocessing %s (%s)...\n',
            (', ').join([str(file) for file in self.input_output['local_paths']]),
            self.model.upper())

        # read in RAP grid from companion netcdf file to deal
        # with cgrib's poor handling of the RAP outer domain grid
        if self.model == 'rap':
            rap_lat_lon_grid = utils.get_rap_lat_lon_grid(Path(self.model_config['lat_lon_grid']))

        # declare empty dict to hold each vertical coordinate's preprocessed data
        preprocessed_data = {}

        # as prescribed in variable_mapping.json: loop through vertical coordinates
        # and the variables we want from each, extracting each variable from GRIB file
        for coord_type,coord_vars in self.model_var_map['pre-process'].items():
            # extract variables that share a vertical coordinate and merge into
            # single Xarray Dataset sharing this vertical coordinate
            var_attrs, merged_vars_by_coord_ds = self._extract_vars_for_vertical_coord(
                coord_type, coord_vars)
            # verify that vertical coordinate has variables
            if len(merged_vars_by_coord_ds.data_vars) == 0:
                logging.warning('No %s variables found! Removing from post-process...\n',
                coord_type)
                del self.model_var_map['post-process'][coord_type]
            else:
                # do any preprocessing that makes sense to do on a variable by variable basis#
                for var_name in merged_vars_by_coord_ds:
                    if self.model == 'rap':
                        # add RAP grid coordinates/dimensions to each variable
                        merged_vars_by_coord_ds[var_name] = utils.add_rap_lat_lon_grid(
                            rap_lat_lon_grid, coord_type, merged_vars_by_coord_ds[var_name])
                    # add variable attributes like units, long names, etc. to each variable
                    merged_vars_by_coord_ds[var_name].attrs = var_attrs[var_name]
                # prep latitude/longitudes as needed, including:
                #    if regular latitude/longitude grid (ie 1D), convert to 2D arrays
                #    recast longitudes from 0 to 360 to -180 to 180 (if needed)
                merged_vars_by_coord_ds = utils.prepare_latitude_longitude(merged_vars_by_coord_ds)
                # calculate rotation grid and add to coordinates
                merged_vars_by_coord_ds = utils.add_rot_grid(merged_vars_by_coord_ds, self.model)
                # add vertical coordinate type to Dataset as attribute:
                merged_vars_by_coord_ds.attrs['vertical_coordinate'] = coord_type
                # for model data on pressure levels (eg GFS), flip the coordinate (so greatest
                # pressure to least) and add pressure variable
                if coord_type == 'isobaricInhPa':
                    merged_vars_by_coord_ds = utils.prep_isobaric_data(merged_vars_by_coord_ds)
                # add merged Xarray Dataset to dictionary with key as coordinate type
                preprocessed_data[coord_type] = merged_vars_by_coord_ds

        return preprocessed_data

    def post_process(self, closest_gridpoint_data: dict[str, xr.Dataset]) -> dict[str, xr.Dataset]:
        """Post-process closest gridpoint model data into a state ready for interpolation

        This includes:
        - Derive any additional fields (model-dependent)
        - Recast winds from grid-relative to earth-relative
        - Scale geopential height on vertical levels so lowest level is 5m, all others relative
        - Perform all unit conversions
        - Interpolate all 3D variables to standard height levels as prescribed by
          new_heights_m in column_extractor.cfg

        :param closest_gridpoint_data: dictionary containing Xarray Datasets of closest data
        :type closest_gridpoint_data: dict
        :returns: Vertical coord-wise dictionary containing Xarray Datasets of closest data
        :rtype: dict
        """

        logging.info('Post-processing the data...\n')

        logging.info('Deriving additional fields from preprocessed data:')

        post_processed_data = closest_gridpoint_data

        # store vertical coordinate name for referencing all 3D work
        vertical_coord = self.model_config['vertical_coord']

        # capture first level of geopotential height, which we consider to be
        # elevation above mean sea level
        # note: this method is physically meaningless for isobaric coordinates
        # so we don't include an elevation field when processing isobaric data (hence the if)
        if vertical_coord != 'isobaricInhPa':
            post_processed_data['surface']['elev'] = calc.get_elevation_above_msl(
                post_processed_data[vertical_coord]['gh'])

        # recast geopotential height on vertical levels to 5m
        post_processed_data[vertical_coord]['gh']  = calc.scale_geopotential_height(
            post_processed_data[vertical_coord]['gh'])

        # calculating mixing ratio at 2m
        post_processed_data['heightAboveGround']['r']  = calc.get_mixing_ratio(
            post_processed_data['heightAboveGround']['q'])

        # calculate mixing ratio on vertical levels from specific humidity
        post_processed_data[vertical_coord]['r'] = calc.get_mixing_ratio(
            post_processed_data[vertical_coord]['q'])

        # calculate virtual temperature on vertical levels from temperature / specific humidity
        post_processed_data[vertical_coord]['tv']  = calc.get_virtual_temp(
            post_processed_data[vertical_coord]['t'],
            post_processed_data[vertical_coord]['r'])

        # calculate vertical velocity on vertical levels from omega, pressure, temperature
        post_processed_data[vertical_coord]['w']  = calc.get_vertical_velocity(
            post_processed_data[vertical_coord]['w'],
            post_processed_data[vertical_coord]['pres'],
            post_processed_data[vertical_coord]['t'])

        # recast u, v components on vertical levels from grid-relative to earth-relative
        # also return wind speed calculated from earth-relative winds
        (post_processed_data[vertical_coord]['u'],
            post_processed_data[vertical_coord]['v'],
            post_processed_data[vertical_coord]['wspd']) = calc.earth_relative_winds(
                post_processed_data[vertical_coord]['u'],
                post_processed_data[vertical_coord]['v'],
                post_processed_data[vertical_coord]['rotation'])

        # recast u, v components at 2m, 10m, 80m from grid-relative to earth-relative
        # (we don't need wind speed at these levels)
        (post_processed_data['heightAboveGround']['u'],
            post_processed_data['heightAboveGround']['v'], *_) = calc.earth_relative_winds(
                post_processed_data['heightAboveGround']['u'],
                post_processed_data['heightAboveGround']['v'],
                post_processed_data['heightAboveGround']['rotation'],
                calc_wind_speed=False)

        logging.info('All additional fields derived. Moving on...\n')

        # interpolate vertical data onto Dave's standardized height levels,
        # as prescribed by new_heights_m in column_extractor.cfg
        # (using gridpoint geopotential height)

        logging.info('Interpolating hybrid level variables onto standard vertical levels:')

        # store geopotential heights for use in interpolation
        geopotential_height = post_processed_data[vertical_coord]['gh']
        # drop geopotential height, because we don't need it
        post_processed_data[vertical_coord] = post_processed_data[vertical_coord].drop_vars(['gh'])

        # parse new heights out of new_heights_m in column_extractor.cfg
        new_heights = list(map(float, self.model_config['new_heights_m'].split(',')))

        # create vertical coordinate Dataset with new heights
        # and drop previous vertical coordinate in the process
        vertically_interpolated_data = post_processed_data[vertical_coord].assign_coords(
            {'height':new_heights})[['height','profiling_site','y','x','distance']]

        # loop through each vertical coordinate variable and perform interpolation
        for var_name,var_data in post_processed_data[vertical_coord].data_vars.items():
            # only interpolate variables that actually have a vertical coordinate
            # (ie not profiling site to model gridpoint distance)
            if vertical_coord in var_data.dims:
                logging.info('--- Interpolating %s', var_name)
                # drop levels with NaNs from variable and geopotential height
                # (this is needed because some GFS variables are available at higher
                # vertical resolution than others
                var_data_for_interp = var_data.dropna(vertical_coord,how='all')
                geopotential_height_for_interp = geopotential_height.sel(
                    {vertical_coord:var_data_for_interp[vertical_coord]}).data
                # this is a Numba-wrapped function that calls numpy's interp()
                # for each gridpoint
                vertically_interpolated_variable = calc.vertical_interpolation(
                    geopotential_height_for_interp,
                    np.array(new_heights),
                    var_data_for_interp.data,
                    np.zeros(tuple([len(new_heights)]) + geopotential_height_for_interp.shape[1:]))
                # store the interpolated variable in the vertical coordinate Dataset
                vertically_interpolated_data[var_name] = xr.DataArray(
                    vertically_interpolated_variable,
                    coords=vertically_interpolated_data.coords,
                    dims=vertically_interpolated_data.dims,
                    attrs=var_data.attrs)
                vertically_interpolated_data[var_name].attrs['vertical_coordinate'] = 'height'

        # replace the vertical level post processed Dataset with the vertically interpolated Dataset
        post_processed_data['height'] = vertically_interpolated_data
        post_processed_data['height'].attrs['vertical_coordinate'] = 'height'

        # get rid of the vertical level post processed Dataset because it's no longer needed
        del post_processed_data[vertical_coord]

        logging.info('All vertical interpolation completed. Moving on...\n')

        logging.info('Performing any necessary unit conversions:')

        # perform any unit conversions from native units to units for output, as prescribed by
        # 'convert_to' field in variable_mapping.json. The native units are extracted from
        # the GRIB file when it's decoded, per ecCodes / local definitions
        for coord_type,coord_vars in self.model_var_map['post-process'].items():
            for var_info in coord_vars:
                if 'convert_to' in var_info.keys():
                    # pass variable to convert to convert_units()
                    post_processed_data[coord_type][var_info['var_name']] = calc.convert_units(
                        post_processed_data[coord_type][var_info['var_name']],
                        var_info['convert_to'])

        logging.info('All unit conversions complete. Moving on...\n')

        return post_processed_data

    def _get_closest_gridpoints(self, preprocessed_data: dict[str, xr.Dataset]) -> dict[str, dict]:
        """Retrieve closest model gridpoint for each profiling site

        :param preprocessed_data: dictionary containing Xarray Datasets of pre-processed data
        :type preprocessed_data: dict
        :returns: dict containing closest y/x indices and distance from for each profiling site
        :rtype: dict
        """

        logging.info('Retrieving closest model gridpoint to each profiling site')

        # store model grid latitude/longitude for non-Xarray calculations
        grid_latitude = preprocessed_data[self.model_config['vertical_coord']]['latitude'].data
        grid_longitude = preprocessed_data[self.model_config['vertical_coord']]['longitude'].data

        # store number of neighbors to extract indices for:
        neighbors = self.model_config.getint('neighbors')

        # create dict to add closest gridpoint + distance between for each profiling site
        closest_gridpoints = {}

        # loop through each profiling site, identify closest model grid point
        for site, site_coords in self.model_profiling_sites.items():
            closest_gridpoints[site] = calc.find_closest_point(
                grid_latitude,
                grid_longitude,
                site_coords)
            # create y index range from closest y index +/- number of neighbors we want to include
            closest_gridpoints[site]['y_range'] = np.arange(
                closest_gridpoints[site]['closest_y'] - neighbors,
                closest_gridpoints[site]['closest_y'] + neighbors + 1)
            # create x index range from closest x index +/- number of neighbors we want to include
            closest_gridpoints[site]['x_range'] = np.arange(
                closest_gridpoints[site]['closest_x'] - neighbors,
                closest_gridpoints[site]['closest_x'] + neighbors + 1)
            # ensure that closest model gridpoint and its neighbors are fully contained
            # within model grid. If not, drop the profiling site
            if (np.logical_or(closest_gridpoints[site]['y_range'].min() < 0,
                closest_gridpoints[site]['x_range'].min() < 0) or
                np.logical_or(closest_gridpoints[site]['y_range'].max() > grid_latitude.shape[0],
                closest_gridpoints[site]['x_range'].max() > grid_latitude.shape[1])):
                logging.info(('--- %sx%s nearest neighborhood to %s is not completely'
                    ' within model grid. Ignoring this site.'), neighbors, neighbors, site)
                del closest_gridpoints[site]

        return closest_gridpoints

    def extract_closest_gridpoint_data(
        self,
        preprocessed_data: dict[str, xr.Dataset]) -> dict[str, xr.Dataset]:
        """Extract model data at each model gridpoint deemed closest to profiling site

        :param preprocessed_data: dictionary containing Xarray Datasets of pre-processed data
        :type preprocessed_data: dict
        :returns: vertical coord-wise dict containing Xarray Datasets for closest model gridpoints
        :rtype: dict
        """

        # retrieve closest gridpoint for each profiling site
        closest_gridpoints = self._get_closest_gridpoints(preprocessed_data)

        # get closest gridpoint data from closest_gridpoints dict
        profiling_sites, distance, closest_y, closest_x, y_range, x_range = map(np.array,
            zip(*[[site] + list(info.values()) for site,info in closest_gridpoints.items()]))

        logging.info('Extracting model data at each closest model gridpoint')

        # create Dataset that contains all relevant closest model gridpoint coordinate information
        # to use when extracting the model data from the original domain
        closest_coords_ds = xr.Dataset(
            data_vars={
                'closest_y':(
                    ['profiling_site'],
                    closest_y,
                    {'units':'unitless','long_name':'Y index of closest gridpoint'}),
                'closest_y_range':(
                    ['profiling_site','y'],
                    y_range,
                    {'units':'unitless','long_name':'Nearest Y neighbors to closest gridpoint'}),
                'closest_x':(
                    ['profiling_site'],
                    closest_x,
                    {'units':'unitless','long_name':'X index of closest gridpoint'}),
                'closest_x_range':(
                    ['profiling_site','x'],
                    x_range,
                    {'units':'unitless','long_name':'Nearest X neighbors to closest gridpoint'}),
                'distance':(
                    ['profiling_site'],
                    distance,
                    {'units':'km','long_name':'Profiling site to closest gridpoint distance'})},
            coords={
                'profiling_site':(['profiling_site'], profiling_sites),
                'y':(['y'], np.arange(y_range.shape[-1])),
                'x':(['x'], np.arange(x_range.shape[-1]))})

        # declare dictionary to add closest gridpoint data to for each vertical coordinate
        closest_gridpoint_data = {}

        # loop through each vertical coordinate and extract data for closest model gridpoints
        for coord_type,coord_data in preprocessed_data.items():
            # extract data for closest model gridpoint +/- neighbors (using their indices)
            closest_gridpoint_data[coord_type] = coord_data.isel(
                y=closest_coords_ds['closest_y_range'],
                x=closest_coords_ds['closest_x_range'])
            # add closest x/y indices and great circle distance as additional coordinate
            closest_gridpoint_data[coord_type] = closest_gridpoint_data[coord_type].assign_coords(
                {
                    'closest_y':closest_coords_ds['closest_y'],
                    'closest_x':closest_coords_ds['closest_y'],
                    'distance':closest_coords_ds['distance']})

        logging.info('Closest model gridpoint data extraction complete. Moving on...\n')

        return closest_gridpoint_data

    def _prepare_interpolated_data_dict(
        self,
        post_processed_data: dict[str, xr.Dataset]) -> dict[str, xr.Dataset]:
        """Prepare the dictionary that will hold vertical coordinate-wise interpolated data

        :param post_processed_data: dictionary containing Xarray Datasets of post-processed_data
        :type post_processed_data: dict
        :returns: vertical coord-wise dict containing empty Xarray datasets by profile site
        :rtype: dict
        """

        # retrieve arrays of profiling site latitude/longitudes
        # (only if they made it into the post-processed data,
        # indicating that they're fully within the model domain)
        site_lats, site_lons = map(
            np.array,
            zip(*[coords for site,coords in
                self.model_profiling_sites.items()
                if site in post_processed_data['height']['profiling_site']]))

        # declare empty dictionary to fill
        interp_data_struct = {}

        # for each vertical coordinate, create an empty Dataset with the correct
        # coordinates and dimensions (profiling site-wise)
        for coord_type,coord_data in post_processed_data.items():
            # instantiate empty interpolated data Dataset, drop y/x dimensions
            interp_data_struct[coord_type] = coord_data.drop_dims(['y','x'])
            for lat_lon_str,lat_lons in zip(('latitude','longitude'), (site_lats,site_lons)):
                # add profiling site latitude/longitudes as coordinates
                interp_data_struct[coord_type] = interp_data_struct[coord_type].assign_coords(
                    {lat_lon_str:xr.DataArray(
                        lat_lons,
                        coords=coord_data['profiling_site'].coords,
                        dims=['profiling_site'],
                        attrs={'units':coord_data[lat_lon_str].attrs['units'],
                            'long_name':f'profiling site {lat_lon_str}'})})
            # remove some coordinates / variables we don't need anymore
            interp_data_struct[coord_type] = interp_data_struct[coord_type].drop_vars(
                ['closest_x', 'closest_y', 'distance'])

        return interp_data_struct

    def interpolate_to_profiling_sites(
        self,
        post_processed_data: dict[str, xr.Dataset]) -> dict[str, xr.Dataset]:
        """Horizontally interpolate post-processed model data to profiling sites of interest

        This includes:
        - Prepare the interpolated data dictionary to be populated with Xarray Datasets
        - Perform 2D interpolation in the horizontal for each vertical coordinate, variable,
          vertical level, and at each profiling site

        :param post_processed_data: dictionary containing Xarray Datasets of post-processed data
        :type post_processed_data: dict
        :returns: dict containing vertical coord-wise Xarray Datasets (interpolated)
        :rtype: dict
        """

        logging.info('Performing horizontal interpolation to profiling sites:')

        # prep the interpolated data Xarray Dataset structure
        interpolated_data = self._prepare_interpolated_data_dict(post_processed_data)

        # store lat/lons of nearest neighbors for interpolation:
        neighbor_lats = post_processed_data['height']['latitude'].data
        neighbor_lons = post_processed_data['height']['longitude'].data

        # store lat/lons of profiling sites for interpolation:
        site_lats = interpolated_data['height']['latitude'].data
        site_lons = interpolated_data['height']['longitude'].data

        # loop through each vertical coordinate, each variable and perform interpolation
        # for each set of profiling sites and their nearest neighbors
        for coord_type,coord_data in post_processed_data.items():
            logging.info('--- Interpolating vertical coordinate type: %s', coord_type)
            # loop through each data variable for this vertical coordinate
            for var_name,var_data in coord_data.data_vars.items():
                logging.info('    Interpolating %s', var_name)
                # perform horizontal interpolation, create Xarray DataArray from the
                # numpy array output. and add to interpolated data Dataset
                interpolated_data[coord_type][var_name] = calc.horizontal_interpolation(
                    neighbor_lats,
                    neighbor_lons,
                    site_lats,
                    site_lons,
                    var_data).assign_coords(interpolated_data[coord_type].coords)

        logging.info('Profiling site interpolation complete. Moving on...\n')

        return interpolated_data

    def write_output_to_netcdf(
        self,
        post_processed_data: dict[str, xr.Dataset],
        interpolated_data: dict[str, xr.Dataset]) -> Path:
        """Write closest model gridpoint data / interpolated data to netCDF

        This includes:
        - Writing temporary CDL file from template CDL file (fill in the placeholders)
        - Creating netCDF file from temporary CDL file (and returning Dataset object to populate)
        - Populating the Dataset object with our data (different handling for different variables)

        Note that we write all netCDF output files to /data/output/... by default. To access this
        file elsewhere, one can mount a local directory into the container as /data/output or write
        the output file to s3 by providing a valid s3 URI to the -o command line argument.

        :param post_processed_data: dict containing vertical coord-wise Xarray Datasets (closest)
        :type post_processed_data: dict
        :param interpolated_data: dict containing vertical coord-wise Xarray Datasets (interped)
        :type interpolated_data: dict
        :returns: filename of local output file
        :rtype: Path
        """

        logging.info('Writing closest model gridpoint / interpolated data to netCDF:')

        # prepare the output netCDF file
        # - create temporary CDL file from template CDL file
        # - write empty netcdf from temporary CDL file, return netCDF object to populate in a bit
        output_dataset, local_filename = file_io.prep_output_file(
            vars(self), interpolated_data, cleanup_cdl=True)

        # note that, in order to retain the original output variable naming conventions, we
        # remap the current variable name to its name in the output file via the `netcdf_mapping`
        # entry for each variable in the `post-process` section for each model in
        # variable_mapping.json

        # write static variables to netCDF file (require slightly different handling):
        # - latitude/longitudes for closest model gridpoints / profiling sites
        # - profiling site names
        # - distance between closest model gridpoints and profiling sites
        # - height bins / depth bins

        logging.info('--- Writing static variables')

        # closest model gridpoint latitudes/longitudes, distances, and height bins
        # (these variables are the same for each vertical coordinate, so just
        # use those from `height`, since this is only coord with height bins)
        for nc_var,xr_var in zip(
            ['mlat','mlon','distance','height'],
            ['latitude','longitude','distance','height']):
            # for the 2D variables (latitude/longitude/height), take all data
            if len(post_processed_data['height'][xr_var].shape) == 1:
                output_dataset[nc_var][:] = post_processed_data['height'][xr_var].data
            # for the 3D variables (distance), pick the center of the nearest neighbors
            # (ie the closest model gridpoint)
            else:
                output_dataset[nc_var][:] = post_processed_data['height'][xr_var].isel(
                    y=self.model_config.getint('neighbors'),
                    x=self.model_config.getint('neighbors')).data

        # profiling site latitudes/longitudes, site names, and depth bins
        # (these variables are the same for each vertical coordinate, so just
        # use those from `depthBelowLandLayer`, since this is only coord with depth bins)
        for nc_var,xr_var in zip(
            ['dlat','dlon','dsite','depth'],
            ['latitude','longitude','profiling_site','depthBelowLandLayer']):
            output_dataset[nc_var][:] = interpolated_data['depthBelowLandLayer'][xr_var].data

        # separate handling for model initialization, because timestamp
        output_dataset['model_initialization'][:] = (
            post_processed_data['height']['time'].data.astype('datetime64[s]').item().timestamp())

        # separate handling for forecast hour / time step, because need to do some math
       # output_dataset['forecast'][:] = calc.forecast_hour_from_timestep(post_processed_data['height']['step'])
        output_dataset['forecast'][:]=3.5 # This will be re-prossessed in bash script using CDO, Xia Sep 2024
        logging.info('--- Writing N-D meteorological variables')

        # write multidimensional meteorological variables to netCDF file:

        # separate handling for elevation above MSL
        # (depending on vertical coordinate type, elevation field might not exist)
        if 'elev' in post_processed_data['surface']:
            output_dataset['elev'][:] = post_processed_data['surface']['elev'].isel(
                surface=0,
                y=self.model_config.getint('neighbors'),
                x=self.model_config.getint('neighbors')).data

        # loop through vertical coordinates / variables in post-process section of
        # model variable map, write to netCDF based on netCDF variable name prescribed
        # by each variable's `netcdf_mapping`
        for coord_type,coord_vars in self.model_var_map['post-process'].items():
            # prep closest model gridpoint data for output netCDF format by selecting
            # nearest neighbor indices and transpose to [profiling site, coordinate]
            # for all variables
            closest_vars = post_processed_data[coord_type].isel(
                y=self.model_config.getint('neighbors'),
                x=self.model_config.getint('neighbors')).transpose(
                    'profiling_site',coord_type)
            # prep interpolated data for output netCDF format by
            # transposing to [profiling site, coordinate]
            interpolated_vars = interpolated_data[coord_type].transpose(
                'profiling_site',
                coord_type)
            # loop through variables for each vertical coordinate
            for var_info in coord_vars:
                # some vertical coordinates have variables with multiple vertical levels that
                # we split into individual variables for the output files
                # (think u at 10m, u at 80m, both in heightAboveGroundLayer Dataset)
                # these are represented by lists in `netcdf_mapping`
                # (ie heightAboveGroundLayer has 2m, 10m, and 80m levels for u
                # which map to netcdf variable names `["N/A", "usfc", "u80m"]` (because no 2m u))
                # if this variable only maps to a single netCDF variable then just write it
                if len(var_info['netcdf_mapping']) == 1:
                    # write output for closest model gridpoint
                    # (netcdf_mapping with 0 appended, per original naming convention)
                    output_dataset[var_info['netcdf_mapping'][0]+'0'][:] = np.expand_dims(
                        closest_vars[var_info['var_name']].data,0)
                    # write output for interpolated variable
                    # (netcdf_mapping with 1 appended, per original naming convention)
                    output_dataset[var_info['netcdf_mapping'][0]+'1'][:] = np.expand_dims(
                        interpolated_vars[var_info['var_name']].data,0)
                # if this variable's levels map to multiple indepedent netCDF variables
                # (it splits into multiple) then loop through each netCDF variable to write
                # (using an integer count to loop through the variable's levels
                else:
                    for nc_var_level,nc_var_name in enumerate(var_info['netcdf_mapping']):
                        # only write variables not listed as `N/A` (meaning that data exists)
                        if nc_var_name != 'N/A':
                            # write output for closest model gridpoint
                            # (netcdf_mapping with 0 appended, per original naming convention)
                            output_dataset[nc_var_name+'0'][:] = np.expand_dims(
                                closest_vars[var_info['var_name']][:,nc_var_level].data,0)
                            # write output for interpolated variable
                            # (netcdf_mapping with 1 appended, per original naming convention)
                            output_dataset[nc_var_name+'1'][:] = np.expand_dims(
                                interpolated_vars[var_info['var_name']][:,nc_var_level].data,0)

        logging.info('Output writing complete. Moving on...\n')

        # close netCDF file
        output_dataset.close()

        return local_filename
