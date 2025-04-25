"""
Description
-----------

Methods for calculating derived fields, rotation coordinates, performing horizontal/vertical
interpolation, etc. This module is imported by, and used throughout, model_data.py.

"""

import logging
from typing import Optional
import pint
import numpy as np
import xarray as xr
from scipy.interpolate import griddata # type: ignore
from numba import jit # type: ignore

# instantiate unit registry
UNIT_REG = pint.UnitRegistry()
# define some extra units we need:
# geopotential meters (which are ~meters in the troposphere)
UNIT_REG.define('geopotential_meter = 1 * meter = gpm')

#### DERIVATIONS, UNIT CONVERSIONS, AND GRID TRANSFORMATION

def get_hrrr_rotation_grid(longitude: xr.DataArray) -> xr.DataArray:
    """Derive HRRR's rotation grid, as per: https://rapidrefresh.noaa.gov/faq/HRRR.faq.html

    We need to convert HRRR winds from grid-relative to earth-relative. That equation requires
    an array that defines the grid rotation (in radians) at each latitude/longitude coordinate
    for a Lambert Conformal projection. This subroutine derives that rotation array, as per:
    https://rapidrefresh.noaa.gov/faq/HRRR.faq.html.

    Note: the rotation equation below has its sign flipped (relative to what's listed at the link).
    This is done to account for the direction of rotation, while still letting us use a single
    set of equations for converting grid-relative winds to earth-relative winds. See
    earth_relative_winds() docString for more information on this.

    :param longitude: DataArray containing HRRR longitude grid
    :type var_: Xarray.DataArray
    :returns: new DataArray containing rotation grid defining grid rotation (in radians)
    :rtype: Xarray.DataArray
    """

    ## define constants, as per: https://rapidrefresh.noaa.gov/faq/HRRR.faq.html
    rotcon_p = 0.622515 # this is sin of lat_tan_p (= 1 for polar stereo)
    lon_xx_p = -97.5 # meridian aligned with cartesian x-axis (in degrees)

    # perform calculation, as per: https://rapidrefresh.noaa.gov/faq/HRRR.faq.html
    rotation = rotcon_p * (longitude - lon_xx_p) * 0.017453
    # flip the sign, to account for the sign change in the
    # grid-relative -> earth-relative wind equations
    # (for rot < 90deg, -sin(rot) = sin(-rot) and cos(rot) = cos(-rot))
    rotation = -1 * rotation

    return rotation.rename('rotation').assign_attrs({
        'standard_name':'rotation',
        'long_name': ("HRRR rotation grid as calculated from equation"
        " defined at https://rapidrefresh.noaa.gov/faq/HRRR.faq.html")})

def get_wfip3_rotation_grid(longitude: xr.DataArray) -> xr.DataArray:
    """Derive WRF WFIP3's rotation grid, adpated from get_hrrr_rotation_grid

    :param longitude: DataArray containing HRRR longitude grid
    :type var_: Xarray.DataArray
    :returns: new DataArray containing rotation grid defining grid rotation (in radians)
    :rtype: Xarray.DataArray
    """
    ## define constants, adapted from: https://rapidrefresh.noaa.gov/faq/HRRR.faq.html
    rotcon_p = 0.650371 # this is sin of lat_tan_p (40.6) (= 1 for polar stereo)
    lon_xx_p = -70.8 # meridian aligned with cartesian x-axis (in degrees)

    # perform calculation, as per: https://rapidrefresh.noaa.gov/faq/HRRR.faq.html
    rotation = rotcon_p * (longitude - lon_xx_p) * 0.017453
    # flip the sign, to account for the sign change in the
    # grid-relative -> earth-relative wind equations
    # (for rot < 90deg, -sin(rot) = sin(-rot) and cos(rot) = cos(-rot))
    rotation = -1 * rotation

    return rotation.rename('rotation').assign_attrs({
        'standard_name':'rotation',
        'long_name': ("HRRR rotation grid as calculated from equation"
        " defined at https://rapidrefresh.noaa.gov/faq/HRRR.faq.html")})

def convert_units(var_: xr.DataArray, to_units: str) -> xr.DataArray:
    """Convert an Xarray DataArray from its native units to to_units via Pint

    Note that to_units can either be true units (ie like "m" for m -> km) or it can be
    "base units", which denotes that we want to convert from current units to the base
    unit of that unit. This is specifically helpful for converting percentages
    (ie 0-100%) to fractions (ie 0-1.0) for cloud cover variables

    :param var_: DataArray to convert units
    :type var_: Xarray.DataArray
    :param to_units: units to convert the DataArray to
    :type to_units: str
    :returns: new DataArray containing the same data, but converted to new units
    :rtype: Xarray.DataArray
    """

    logging.info('--- Converting %s (%s) from %s to %s',
        var_.name, var_.attrs['standard_name'], var_.units, to_units)

    # make copy to avoid overwriting original
    var_to_convert = var_.copy(deep=False)

    # perform conversion
    # this may be to another unit like m -> km
    # or a base unit like percentage (0-100%) -> fraction (0-1.0)
    if to_units == 'base units':
        converted_var = UNIT_REG.Quantity(var_to_convert.data,var_to_convert.units).to_base_units()
    else:
        converted_var = UNIT_REG.Quantity(var_to_convert.data,var_to_convert.units).to(to_units)

    # replace data attribute of var_to_convert with the magnitude (ie no units)
    # of the converted variable
    var_to_convert.data = converted_var.magnitude

    # return var_to_convert with new units assigned to units attribute
    return var_to_convert.assign_attrs({'units':str(converted_var.units)})

def get_mixing_ratio(specific_humidity: xr.DataArray) -> xr.DataArray:
    """Calculate mixing ratio from specific humidity

    This subroutine leverages the derivation in Hobbs and Wallace (2006) as:
    mixing_ratio = specific_humidity / (1 - specific_humidity). The returned
    mixing_ratio DataArray is dimensionless.

    :param specific_humidity: DataArray containing the specific humidity variable
    :type specific_humidity: Xarray.DataArray
    :returns: new DataArray containing mixing ratio in dimensionless units
    :rtype: Xarray.DataArray
    """

    logging.info('--- Deriving mixing ratio (%s)', specific_humidity.attrs['GRIB_typeOfLevel'])

    # ensure that specific humidity is dimensionless
    specific_humidity = convert_units(specific_humidity, 'dimensionless')

    # perform derivation, using equation in Hobbs and Wallace (2006)
    mixing_ratio = specific_humidity / (1 - specific_humidity)

    # set DataArray attributes (units are already set when returned as Quantity)
    return mixing_ratio.assign_attrs({
        'standard_name':'mixing_ratio',
        'long_name':("Mixing ratio calculated from specific humidity following"
        " the derivation in Hobbs and Wallace (2006)"),
        'vertical_coordinate':specific_humidity.attrs['vertical_coordinate'],
        'units':'dimensionless'})

def get_virtual_temp(temperature: xr.DataArray, mixing_ratio: xr.DataArray) -> xr.DataArray:
    """Calculate virtual temperature from temperature and mixing ratio

    This subroutine uses the derivation defined in Hobbs and Wallace (2006) as:
    T_v = temperature * ((mixing_ratio + epsilon) / (epsilon * (1 + mixing_ratio))),
    where epsilon is the molecular weight ratio between water vapor and dry air.
    The returned virtual_temp DataArray has units of K.

    :param temperature: DataArray containing the temperature variable
    :type temperature: Xarray.DataArray
    :param mixing_ratio: DataArray containing the mixing ratio variable
    :type mixing_ratio: Xarray.DataArray
    :returns: new DataArray containing virtual temperature in units of K
    :rtype: Xarray.DataArray
    """

    logging.info('--- Deriving virtual temperature (%s)', temperature.attrs['GRIB_typeOfLevel'])

    # define molecular weight ratio between water vapor and dry air (epsilon)
    epsilon = 0.6219569100577033

    # ensure that temperature is in K
    temperature = convert_units(temperature, 'K')
    # ensure that mixing ratio is dimensionless
    mixing_ratio = convert_units(mixing_ratio, 'dimensionless')

    # perform derivation, using derivation in Hobbs and Wallace (2006)
    virtual_temp = temperature * ((mixing_ratio + epsilon) / (epsilon * (1 + mixing_ratio)))

    return virtual_temp.assign_attrs({
        'standard_name':'virtual_temperature',
        'long_name':("Virtual temperature calculated from temperature and mixing"
        " ratio using derivation found in Hobbs and Wallace (2006)"),
        'vertical_coordinate':temperature.attrs['vertical_coordinate'],
        'units':'K'})

def _get_density(pressure: xr.DataArray, temperature: xr.DataArray) -> xr.DataArray:
    """Calculate density from pressure and temperature, assuming dry air

    This follows the equation rho = pressure / (Rd * temperature), where rho is density and Rd
    is the dry air constant.

    :param pressure: DataArray containing the pressure variable
    :type pressure: Xarray.DataArray
    :param temperature: DataArray containing the temperature variable
    :type temperature: Xarray.DataArray
    :returns: new DataArray containing density in units __
    :rtype: Xarray.DataArray
    """

    # define dry air constant (Rd) as quotient of molar gas constant to dry air molecular weight
    # in units J / kg**-1 K**-1
    dry_air_constant =  8.314462618 / 28.96546e-3

    return pressure / (dry_air_constant * temperature)

def get_vertical_velocity(
    omega: xr.DataArray,
    pressure: xr.DataArray,
    temperature: xr.DataArray) -> xr.DataArray:
    """Calculate vertical velocity from pressure vertical velocity, pressure, and temperature

    This subroutine follows the derivation outlined in Hobbs and Wallace (2006), as:
    vertical_vel = (omega / (-g * rho)), g is the gravitationl constant,
    rho = pressure / (Rd * temperature), and Rd is the dry air constant. This assumes
    dry air and doesn't use the virtual temperature approximation. The returned vertical_vel
    DataArray has units of m/s.

    :param omega: DataArray containing the omega variable
    :type omega: Xarray.DataArray
    :param pressure: DataArray containing the pressure variable
    :type pressure: Xarray.DataArray
    :param temperature: DataArray containing the temperature variable
    :type temperature: Xarray.DataArray
    :returns: new DataArray containing vertical velocity in units m/s
    :rtype: Xarray.DataArray
    """

    logging.info('--- Deriving vertical velocity (%s)', omega.attrs['GRIB_typeOfLevel'])

    # define Earth's gravitational constant (g) in m / s**2
    gravitational_constant = 9.80665

    # ensure that omega is in Pa/s
    omega = convert_units(omega, 'Pa/s')
    # ensure that pressure is in Pa
    pressure = convert_units(pressure, 'Pa')
    # ensure that temperature is in K
    temperature = convert_units(temperature, 'K')

    # calculate density (rho) from pressure and temperature
    density = _get_density(pressure, temperature)

    # perform derivation, following approximation in Hobbs and Wallace (2006)
    vertical_vel = omega / (-gravitational_constant * density)

    return vertical_vel.assign_attrs({
        'standard_name':'vertical_velocity',
        'long_name':("Vertical velocity calculated from omega, temperature, and"
        " pressure using derivation in Hobbs and Wallace (2006)"),
        'vertical_coordinate':omega.attrs['vertical_coordinate'],
        'units':'m/s'})

def earth_relative_winds(
    u_comp: xr.DataArray,
    v_comp: xr.DataArray,
    rotation: xr.DataArray,
    calc_wind_speed: bool = True) -> tuple[xr.DataArray, xr.DataArray, Optional[xr.DataArray]]:
    """Convert u, v components of the wind from grid-relative to earth-relative

    Native u, v components of the wind are in grid-relative terms. They need converted to
    earth-relative. This is done via the following equations:

    u_earth = cos(rotation) * u_grid - sin(rotation) * v_grid
    v_earth = cos(rotation) * v_grid + sin(rotation) * u_grid

    Where rotation is a 2D array with same size as latitude/longitude and is a function of the
    model grid. The returned u_earth_relative and v_earth_relative DataArrays have units of m/s.
    If calc_wind_speed is True, calculate the wind speed from earth-relative u, v components.
    This is the default behavior.

    Note that these equations are slightly different for HRRR, because the rotation direction
    is different (note the flipped signs):

    u_earth = cos(rotation) * u_grid + sin(rotation) * v_grid
    v_earth = cos(rotation) * v_grid - sin(rotation) * u_grid

    However, we can still use the first set of equations by flipping the sign on the HRRR rotation
    grid in get_hrrr_rotation_grid() because the following identities hold for small rotation
    angles (rotation < 90deg):

    cos(-1 * rotation) = cos(rotation)
    sin(-1 * rotation) = -1 * sin(rotation)

    :param u_comp: DataArray containing the u variable
    :type u_comp: Xarray.DataArray
    :param v_comp: DataArray containing the v variable
    :type v_comp: Xarray.DataArray
    :param rotation: DataArray containing the rotation variable
    :type rotation: Xarray.DataArray
    :param calc_wind_speed: whether to calculate and return the wind speed
    :type calc_wind_speed: bool
    :returns: tuple containing earth-relative u, v components / speed of the wind in units m/s
    :rtype: (Xarray.DataArray, Xarray.DataArray, Xarray.DataArray)
    """

    logging.info('--- Converting u, v components from grid-relative to earth-relative (%s)',
        u_comp.attrs['GRIB_typeOfLevel'])

    # ensure that u component is in m/s
    u_comp = convert_units(u_comp, 'm/s')
    # ensure that v component is in m/s
    v_comp = convert_units(v_comp, 'm/s')

    # perform derivation for u component of the wind
    u_earth_relative =  u_comp * np.cos(rotation) - v_comp * np.sin(rotation)

    # set DataArray attributes (units are already set when returned as Quantity)
    u_earth_relative = u_earth_relative.assign_attrs({
        'standard_name':'u_earth_relative',
        'long_name':("u component of the wind in earth-relative coordinates,"
        " as calculated via sin(rotation) * v_grid + cos(rotation) * u_grid"),
        'vertical_coordinate':u_comp.attrs['vertical_coordinate'],
        'units':'m/s'})

    # perform derivation for v component of the wind
    v_earth_relative = v_comp * np.cos(rotation) + u_comp * np.sin(rotation)

    # set DataArray attributes (units are already set when returned as Quantity)
    v_earth_relative = v_earth_relative.assign_attrs({
        'standard_name':'v_earth_relative',
        'long_name':("v component of the wind in earth-relative coordinates,"
        " as calculated via cos(rotation) * v_grid - sin(rotation) * u_grid"),
        'vertical_coordinate':v_comp.attrs['vertical_coordinate'],
        'units':'m/s'})

    # if you asked for the wind speeds, calculate them and return:
    if calc_wind_speed:

        logging.info('--- Calculating wind speed from u, v components (%s)',
            u_comp.attrs['GRIB_typeOfLevel'])

        # this isn't technically earth-relative, since wind speed is conserved regardless
        # of coordinate system (but we'll name it so for conventions' sake)
        wind_speed_relative = np.hypot(u_earth_relative, v_earth_relative)

        # set DataArray attributes (units are already set when returned as Quantity)
        wind_speed_relative = wind_speed_relative.assign_attrs({        # type: ignore
            'standard_name':'wind_speed_earth_relative',                # type: ignore
            'long_name':"wind speed using Numpy's hypot()",             # type: ignore
            'vertical_coordinate':u_comp.attrs['vertical_coordinate'],  # type: ignore
            'units':'m/s'})                                             # type: ignore

        return u_earth_relative, v_earth_relative, wind_speed_relative # type: ignore

    return u_earth_relative, v_earth_relative # type: ignore

def get_elevation_above_msl(geopotential_height: xr.DataArray) -> xr.DataArray:
    """Retrieve elevation above mean sea level from geopotential height

    Per convention, the first level of geopotential heights are considered an approximation
    for elevation above mean sea level.

    :param geopotential_height: DataArray containing the geopotential height variable
    :type geopotential_height: Xarray.DataArray
    :returns: new DataArray containing elevation above mean sea level
    :rtype: Xarray.DataArray
    """

    logging.info('--- Extracting elevation above MSL from geopotential height (first level)')

    # extract first level of geopotential height, which we consider surface elevation above MSL
    elevation = geopotential_height[0,:,:,:]

    # get rid of geopotential height attributes and add new attributes for elevation
    elevation.attrs = {
        'standard_name':'elevation',
        'long_name':("elevation above mean sea level adjusted for gravity"
        " (first level of geopotential height)"),
        'vertical_coordinate':'surface',
        'units':'m'}

    return elevation.expand_dims('surface')

def scale_geopotential_height(geopotential_height: xr.DataArray) -> xr.DataArray:
    """Scale geopotential height field such that first level is at 5m and all other levels relative

    Geopotential heights on vertical coordinates are scaled such that gh_n = (gh_n - gh_0) + 5. The
    geopotential height field is what's eventually used as the vertical coordinate to interpolate
    against, relative to the prescribed new levels. This is done specifically to stay in line with
    the original workflow.

    :param geopotential_height: DataArray containing the geopotential height variable
    :type geopotential_height: Xarray.DataArray
    :returns: new DataArray containing geopotential height scaled as prescribed in units m
    :rtype: Xarray.DataArray
    """

    logging.info('--- Scaling geopotential heights relative to 5m')

    # store vertical coordinate to add to new scaled variable attributes
    vertical_coordinate = geopotential_height.attrs['vertical_coordinate']

    # ensure that geopotential heights are in m
    geopotential_height = convert_units(geopotential_height, 'm')

    # scale geopotential height relative to level_0 = 5m
    geopotential_height = (geopotential_height -
        geopotential_height.isel({vertical_coordinate:0})) + 5

    return geopotential_height.assign_attrs({
        'standard_name':'geopotential_height',
        'long_name':('Geopotential height recasted such that the lowest'
        ' model level is 5m and all other vertical levels are relative to that'),
        'vertical_coordinate':vertical_coordinate,
        'units':'m'})

#### GEOGRAPHIC MATH (GREAT CIRCLE DISTANCES, CLOSEST POINTS, ETC.)

@jit(nopython=True)
def _haversine(
    grid_latitude: np.ndarray,
    grid_longitude: np.ndarray,
    site_latitude: float,
    site_longitude: float) -> np.ndarray:
    """
    Calculate the great circle distance in kilometers between two points on the earth

    This method uses the Haversine formula. For more information, see:
    https://en.wikipedia.org/wiki/Haversine_formula
    We chose Haversine over the more accurate Vincenty method (~0.5%) because it's more than
    an order of magnitude faster and yields the same results for closest model gridpoints. This
    method is juiced with Numba!

    :param grid_latitude: array containing model latitude grid
    :type grid_latitude: numpy.ndarray
    :param grid_longitude: array containing model longitude grid
    :type grid_longitude: numpy.ndarray
    :param site_latitude: latitude for profiling site of interest
    :type site_latitude: float
    :param site_longitude: longitude for profiling site of interest
    :type site_longitude: float
    :returns: array containing great circle distance between profiling site and model gridpoints
    :rtype: numpy.ndarray
    """

    # mean radius of the earth (km) - per First Course in Atmospheric Radiation, 2nd Edition
    earth_radius = 6373.

    # haversine formula
    delta_longitude = np.radians(site_longitude) - np.radians(grid_longitude)
    delta_latitude = np.radians(site_latitude) - np.radians(grid_latitude)
    great_circle_dist = (np.sin(delta_latitude/2.)**2
        + np.cos(np.radians(grid_latitude))
        * np.cos(np.radians(site_latitude))
        * np.sin(delta_longitude/2.)**2)
    central_angle = 2 * np.arcsin(np.sqrt(great_circle_dist))

    # return distance between profiling site of interest and each model gridpoint (km)
    return central_angle * earth_radius

def find_closest_point(
    grid_latitude: np.ndarray,
    grid_longitude: np.ndarray,
    site_coords: tuple[float, float]) -> dict:
    """Identify closest model grid point to profiling site via great circle distance

    Method performs great circle calculations to compute the distance between profiling site and
    each lat/lon point on model grid using the Haversine formula. From this, we deduce the closest
    model gridpoint to the profiling site and return the distance between / indices for this grid
    point.

    :param grid_latitude: array containing model latitude grid
    :type grid_latitude: numpy.ndarray
    :param grid_longitude: array containing model longitude grid
    :type grid_longitude: numpy.ndarray
    :param site_coords: tuple containing the (latitude, longitude) for profiling site of interest
    :type site_coords: tuple
    :returns: dict containing the indices / cartesian distance to closest model gridpoint
    :rtype: dict
    """

    # expand site coords into latitude/longitude
    site_latitude, site_longitude = site_coords

    # calculate distance between site lat/lon and each lat/lon point on model grid
    # uses Haversine formula (order of magnitude faster than Vincenty, yields results within 0.5%)
    distances = _haversine(
        grid_latitude,
        grid_longitude,
        np.tile(site_latitude, grid_latitude.shape),
        np.tile(site_longitude, grid_longitude.shape))

    # retrieve minimum distance (ie, closest point) in km
    min_distance = np.min(distances)

    # retrieve indices at minimum distance (ie, closest point)
    closest_indices = np.unravel_index(np.argmin(distances), distances.shape)

    # return dict containing distance and closest y/x indices
    return ({'distance_between':min_distance} |
        dict(zip(('closest_y','closest_x'), closest_indices)))

#### INTERPOLATION: VERTICAL, HORIZONTAL, ETC

@jit(nopython=True)
def vertical_interpolation(
    current_heights: np.ndarray,
    new_heights: np.ndarray,
    data_to_interpolate: np.ndarray,
    output_array: np.ndarray) -> np.ndarray:
    """Perform 1D interpolation in the vertical for each model gridpoint (juiced with Numba!)

    We do this to resample/interpolate the vertical coordinate variable to a standard vertical
    height coordinate, as prescribed by new_heights_m in column_extractor.cfg. For each gridpoint,
    the geopotential height profile is used as the current height coordinate, which we use to
    interpolate to new height coordinate new_heights_m. This function is wrapped with Numba,
    which compiles it as C and enables ~20x increase in performance.

    :param current_heights: array of geopotential height profiles at each model gridpoint
    :type current_heights: numpy.ndarray
    :param new_heights: array of standardized heights, as prescribed in configuration file
    :type new_heights: numpy.ndarray
    :param data_to_interpolate: array containing data for variable we want to interpolate
    :type data_to_interpolate: numpy.ndarray
    :param output_array: null (zeros) array to add interpolated profiles to
    :type output_array: numpy.ndarray
    :returns: array containing previously empty output array populated with interpolated profiles
    :rtype: numpy.ndarray
    """

    # loop through each profiling site
    for site_index in range(current_heights.shape[1]):
        # loop through y coordinates of grid
        for y_index in range(current_heights.shape[2]):
            # loop through x coordinates of grid
            for x_index in range(current_heights.shape[3]):
                # perform 1D interpolation on the vertical profile defined by
                # [vertical, site_index, y_index, x_index]
                output_array[:, site_index, y_index, x_index] = np.interp(
                    new_heights,
                    current_heights[:, site_index, y_index, x_index],
                    data_to_interpolate[:, site_index, y_index, x_index])

    return output_array

def horizontal_interpolation(
    neighbor_lats: np.ndarray,
    neighbor_lons: np.ndarray,
    target_lats: np.ndarray,
    target_lons: np.ndarray,
    variable: xr.DataArray) -> xr.DataArray:
    """Perform 2D interpolation in the horizontal for each vertical level, each profiling site

    This subroutine uses Scipy's griddata() function to perform linear interpolation. For more
    information, see:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.griddata.html.
    (unfortunately not juiced with Numba due to no Scipy integration)

    :param neighbor_lats: array of latitude neighborhood for each closest model gridpoint
    :type neighbor_lats: numpy.ndarray
    :param neighbor_lons: array of longitude neighborhood for each closest model gridpoint
    :type neighbor_lons: numpy.ndarray
    :param target_lats: array of latitudes at each profiling site
    :type target_lats: numpy.ndarray
    :param target_lons: array of longitudes at each profiling site
    :type target_lons: numpy.ndarray
    :param variable: DataArray containing post processed data for variable to be interpolated
    :type variable: xarray.DataArray
    :returns: DataArray containing variable interpolated to each profiling site
    :rtype: xarray.DataArray
    """

    # load variable's data into memory
    variable_to_interpolate = variable.data

    # create array that will be filled with interpolated data
    interpolated_variable = np.full(
        (variable[variable.attrs['vertical_coordinate']].size,variable['profiling_site'].size),
        np.nan)

    # loop through the levels associated with this vertical coordinate
    for level_idx in range(variable[variable.attrs['vertical_coordinate']].size):
        # some levels don't have any data and are full of NaNs (ie 10m/80m RH, q, etc.)
        # so we save some time by skipping these
        if not np.isnan(variable_to_interpolate[level_idx,:,:,:]).all():
            # if we have data, loop through each profiling site and interpolate
            for site_idx in range(variable['profiling_site'].size):
                interpolated_variable[level_idx,site_idx] = griddata(
                    (neighbor_lons[site_idx,:,:].ravel(),
                    neighbor_lats[site_idx,:,:].ravel()),
                    variable_to_interpolate[level_idx,site_idx,:,:].ravel(),
                    (target_lons[site_idx],target_lats[site_idx]))

    # return as xarray DataArray
    return xr.DataArray(
            interpolated_variable,
            dims=[variable.attrs['vertical_coordinate'],'profiling_site'],
            attrs=variable.attrs)

#### MISCELLANEOUS MATH

def forecast_hour_from_timestep(timestep: xr.DataArray) -> float:
    """Calculate float forecast hour (since initialization) from timestep in nanoseconds

    Native form of timestep is in nanoseconds since initialization. So, for example,
    forecast hour 3 is given as 10800000000000 nanoseconds. Here we convert from nanoseconds
    to float hours from initialization time.

    :param timestep: timestep since initialization in nanoseconds
    :type timestep: Xarray.DataArray
    :returns: float forecast hour
    :rtype: float
    """

    return float(timestep.data / np.timedelta64(1, 'h'))
