"""
Description
-----------

Methods for writing extracted columns data to netCDF files. This includes creating temporary
CDL file from template CDL file, instantiating empty netCDF file from temporary CDL file,
writing data to netCDF Dataset object, and cleaning up the temporary CDL. Also included
are methods for interacting with s3 (whether for fetching input file or pushing output file)

"""

import os
import logging
import re
import json
import configparser
from typing import Any, Optional
from pathlib import Path
import xarray as xr
import numpy as np
from netCDF4 import Dataset # type: ignore # pylint: disable=no-name-in-module
import botocore # type: ignore
from boto3 import client # type: ignore

#### HANDLES INITIAL CONFIGURATION/MODEL FILES

def file_exists(file_: Path) -> None:
    """Check if file exists via pathlib.Path.is_file(). If not, raise FileNotFoundError

    :param file_: the file we want to check existence for
    :type file_: Path
    :returns: {None}
    :rtype: {None}
    """

    if not file_.is_file():
        raise FileNotFoundError(f"{file_.name} doesn't exist. Exiting...\n")

def read_json(json_file: Path) -> dict:
    """Read in json using ``json.load`` method , if the json exists (log and exit otherwise)

    :param json_file: the json file path
    :type json_file: Path
    :returns: nested dictionary containing json contents
    :rtype: dict
    """

    # verify that file exists
    file_exists(json_file)

    logging.info('Reading %s...\n', json_file.name)

    # read in json file
    with open(json_file,'r',encoding='utf-8') as json_file_to_read:
        json_contents_dict = json.load(json_file_to_read)

    return json_contents_dict

def read_config(config_file: Path) -> configparser.ConfigParser:
    """Parse configuration file contents, if the config file exists (log and exit otherwise)

    :param config_file: path to configuration file
    :type config_file: Path
    :returns: parsed contents of configuration file
    :rtype: configparser.ConfigParser
    """

    # verify that file exists
    file_exists(config_file)

    logging.info('Reading %s...\n', config_file.name)

    # read in config file
    config = configparser.ConfigParser(allow_no_value=True)
    config.read(config_file)

    return config

def read_grib(local_file_path: Path, filter_keys: dict[str, str]) -> xr.Dataset:
    """Read input file; extract variable on vertical coordinate based on filter keys

    Since GRIB files can contain arbitrary heterogenous combinations of variables / vertical
    coordinates / step types, etc., it's better to extract each variable one at a time to
    ensure that there are no collisions. As such, we treat each variable as a unique
    combination of GRIB short name (ie `t`), vertical coordinate type (ie `hybrid`), and
    various other GRIB keys. This is done via cfgrib's `filter_by_keys` option which is
    passed via xarray.open_dataset's `backend_kwargs` argument.

    :param local_file_path: local path to file to process
    :type local_file_path: pathlib.Path
    :param filter_keys: dictionary containing all variable info we want to use as filter keys
    :type filter_keys: dict
    :returns: Dataset containing the variable we wanted, for vertical coordinate / step type
    :rtype: xarray.Dataset
    """

    logging.info('--- %s',
        ''.join([f"{key+':':<20s}{value:<25s}" for key,value in filter_keys.items()]))

    # read GRIB file for variable specified by filter_keys using cfgrib as backend
    # (ignore any .idx files that may / may not exist)
    return xr.open_dataset(
        local_file_path,
        engine='cfgrib',
        backend_kwargs={'filter_by_keys':filter_keys})

#### S3 COMPONENTS FOR FETCHING SOURCE FILES AND PUTTING OUTPUT FILES

def validate_s3_uri(s3_uri: str, endpoint_url: Optional[str]) -> Optional[re.Match]:
    """Verify that S3 URI matches format: s3://<bucket>/<key 1>/<key 2>/...

    :param s3_uri: S3 path to verify
    :type s3_uri: str
    :param endpoint_url: S3 endpoint URL or None
    :type endpoint_url: str or None
    :returns: re.Match object containing bucket and key structure or None
    :rtype: re.Match or None
    """

    # check s3_uri to see if it's a correctly formatted S3 path
    # if it is, parse out bucket/key and return
    if parsed_uri := re.match(r's3:\/\/(?P<BUCKET>.+?)\/(?P<KEY>.+)', s3_uri):
        logging.info('%s is a valid S3 URI', s3_uri)
        logging.info('--- %-20s: %s', 'Bucket is', parsed_uri.group('BUCKET'))
        logging.info('--- %-20s: %s', 'Key(s) is', parsed_uri.group('KEY'))
        logging.info('--- %-20s: %s\n', 'Endpoint URL is', endpoint_url)
        return parsed_uri

    raise RuntimeError("Specified S3 URI isn't of format s3://<bucket>/<key 1>/<key 2>/...)\n")

def _instantiate_s3_client(endpoint_url: Optional[str]) -> client:
    """Instantiate S3 client via boto3.

    Handles three access cases:
    - Internal /public via Minio endpoint URL and unsigned request from on-prem
    - GSL s3 buckets using task IAM role from inside AWS ECS tasks
    - Public-facing s3 buckets (ie NODD) that don't need credentials
      from inside AWS ECS tasks or on-prem

    :param endpoint_url: optional URL to endpoint (if accessing internal /public via Minio)
    :type endpoint_url: str or None
    :returns: S3 client via boto3
    :rtype: client
    """

    # if endpoint URL is provided, pass it in and do unsigned request
    # endpoint URL used when fetching internal /public file via Minio
    if endpoint_url:
        s3_client = client(
            's3',
            endpoint_url=endpoint_url,
            config=botocore.config.Config(signature_version=botocore.UNSIGNED))
    # if AWS_CONTAINER_CREDENTIALS_RELATIVE_URI not in environment and no endpoint URL, do unsigned
    # likely working with public-facing s3 bucket that doesn't need credentials
    elif "AWS_CONTAINER_CREDENTIALS_RELATIVE_URI" not in os.environ:
        s3_client = client(
            's3',
            config=botocore.config.Config(signature_version=botocore.UNSIGNED))
    # if AWS_CONTAINER_CREDENTIALS_RELATIVE_URI in environment, then running inside AWS ECS task
    # boto3 automatically handles task IAM role credentialling
    else:
        s3_client = client('s3')

    return s3_client

def download_from_s3(parsed_uri: re.Match, endpoint_url: Optional[str]) -> Optional[Path]:
    """Download file via S3 protocol using AWS's boto3 package

    :param parsed_uri: parsed S3 URI containing bucket and key information
    :type parsed_uri: re.Match
    :param endpoint_url: optional URL to endpoint (if accessing internal /public via Minio)
    :type endpoint_url: str or None
    :returns: path to local file if download was successful
    :rtype: Path or None
    """

    # downloaded file will have same filename as remote file and will land in /data/input
    local_filename = f"/data/input/{parsed_uri['KEY'].split('/')[-1]}"

    logging.info('Attempting to download file to: %s', local_filename)

    # instantiate S3 client via boto3
    s3_client = _instantiate_s3_client(endpoint_url)

    # download the file at the bucket/key provided via file_location / parsed_uri
    s3_client.download_file(
        Bucket=parsed_uri['BUCKET'],
        Key=parsed_uri['KEY'],
        Filename=local_filename)

    logging.info('--- Download successful\n')

    return Path(local_filename)

def push_output_to_s3(local_filename: Path, s3_uri: str, endpoint_url: Optional[str]) -> None:
    """Upload file via S3 protocol using AWS's boto3 package

    :param local_filename: path to local file to upload to S3
    :type local_filename: Path
    :param s3_uri: path to upload file to on S3 (includes bucket/keys)
    :type s3_uri: str
    :param endpoint_url: optional URL to endpoint (if accessing internal /public via Minio)
    :type endpoint_url: str
    :returns: {None}
    :rtype: {None}
    """

    # only attempt push if S3 URI is correctly formatted
    if parsed_uri := validate_s3_uri(s3_uri, endpoint_url):

        # s3 file will have same filename as local file, but land at bucket / key
        s3_key_path = f"{parsed_uri['KEY']}/{local_filename.name}"

        logging.info('Attempting to upload file to: %s', s3_key_path)

        # instantiate S3 client via boto3
        s3_client = _instantiate_s3_client(endpoint_url)

        # upload the file to the bucket/key provided
        s3_client.upload_file(
            Bucket=parsed_uri['BUCKET'],
            Key=s3_key_path,
            Filename=str(local_filename))

        logging.info('--- Upload successful\n')

#### METHODS FOR WRITING OUTPUT FILES

def _get_output_filename(
    config_vars: dict[str, Any],
    cycle_time: str,
    forecast_hour: str) -> tuple[Path, Path]:
    """Construct temporary CDL / netCDF filenames from output path, cycle time, and forecast hour

    Note that we write all netCDF output files to /data/output/... by default. To access this file
    elsewhere, one can mount a local directory into the container as /data/output or write the
    output file to s3 by providing a valid s3 URI to the -o command line argument.

    :param config_vars: dict containing all attributes from ProcConfig class
    :type config_vars: dict
    :param cycle_time: Cycle time for the file being processed
    :type cycle_time: str
    :param forecast_hour: Forecast hour for the file being processed
    :type forecast_hour: str
    :returns: tuple containing full temporary CDL filename and output netCDF filename
    :rtype: (pathlib.Path, pathlib.Path)
    """

    # construct time string for filenames
    time_str = cycle_time +  forecast_hour.zfill(2)

    # construct output filename root
    output_filename_root = f"profiler_sitesF.ncep_{config_vars['model']}.{time_str}"

    # construct full filenames (one for .nc and one for .cdl)
    # note that we write all output files to /data/output/... by default
    cdl_filename, nc_filename = [
        Path('/data/output') / f'{output_filename_root}.{ext}' for ext in ['cdl','nc']]

    return cdl_filename, nc_filename

def _replace_cdl_content(cdl_template_path: Path, metadata: dict[str, Any]) -> str:
    """Replace metadata placeholders in template CDL with relevant metadata, save to temporary CDL

    For flexibility's sake, we have placeholders for the number of profiling sites, height grid,
    as well as relevant metadata like the file being processed, model initialization, current
    code version tag, etc. There are placeholders for each in the template CDLs. This function
    reads in the template CDL, replaces the placeholders with relevant metadata, and saves
    as temporary CDL used to write empty netCDF.

    :param cdl_template_path: path to CDL template for current model
    :type cdl_template_path: pathlib.Path
    :param metadata: dict containing metadata to insert into temporary CDL
    :type metadata: dict
    :returns: formatted string representing the temporary CDL to be saved to disk
    :rtype: str
    """

    logging.info('--- %-30s: %s', 'Modifying CDL template', cdl_template_path)

    # read in the template CDL file
    with open(cdl_template_path, "r", encoding='utf-8') as cdl_template:
        cdl_template_content = cdl_template.read()

    # store template CDL file contents
    modified_cdl_content = cdl_template_content

    # loop through each placeholder string, and replace it with the current metadata
    for str_to_replace,replace_with in metadata.items():
        if str_to_replace in modified_cdl_content:
            modified_cdl_content = modified_cdl_content.replace(str_to_replace,replace_with)

    return modified_cdl_content

def _write_temporary_cdl(
    config_vars: dict[str, Any],
    interpolated_data: dict[str, xr.Dataset]) -> tuple[Path, Path]:
    """Write the temporary CDL we just populated to disk

    This function constructs the CDL/netCDF filenames, inserts relevant metadata into our
    new temporary CDL, and writes it to disk. This temporary CDL is ultimately what is used
    to write the empty netCDF that will be populated with our processed data.

    :param config_vars: dict containing all attributes from ProcConfig class
    :type config_vars: dict
    :param interpolated_data: dict containing vertical coord-wise Xarray Datasets (interpolated)
    :type interpolated_data: dict
    :returns: tuple containing full temporary CDL filename and output netCDF filename
    :rtype: (pathlib.Path, pathlib.Path)
    """

    # store model initialization and timesteps
    model_init = interpolated_data['height']['time']
    model_step = 1800000000      # for debug, Xia, Sep 9, 2024
    # model_step = interpolated_data['height']['step']

    # handle source file metadata
    if None not in config_vars['input_output']['remote_paths']:
        source_file = (', ').join(config_vars['input_output']['remote_paths'])
    else:
        source_file = (', ').join(
            [str(file) for file in config_vars['input_output']['local_paths']])

    # construct metadata dict to poke into CDL:
    metadata = {
        'N_SITES':str(interpolated_data['height'].sizes['profiling_site']),
        'N_HEIGHT_BINS':str(interpolated_data['height'].sizes['height']),
        'N_DEPTH_BINS':str(interpolated_data['depthBelowLandLayer'].sizes['depthBelowLandLayer']),
        'SOURCE_FILE':f'"{source_file}"',
        'MODEL_INIT':f'"{model_init.dt.strftime("%Y-%m-%dT%H:%M:%S").data}"',
        'CODE_VERSION':f'"{os.environ["CODE_VERSION"]}"',
        'FILE_CREATION_TIME':f'"{config_vars["start_time"].strftime("%Y-%m-%dT%H:%M:%S")}"'}

    # construct filenames for temporary CDL and netCDF files
    cdl_filename, nc_filename = _get_output_filename(
        config_vars,
        str(model_init.dt.strftime('%Y%j.%H').data),
        str('subhour'))   # Str needs to be replaced with the real sub-hourly forecast hour, Xia, Sep 2024

    # insert relevant metadata into contents of CDL template
    modified_cdl_content = _replace_cdl_content(
        config_vars['input_output']['model_output_template_path'],
        metadata)

    logging.info('--- %-30s: %s', 'Writing CDL', cdl_filename)

    # write CDL that contains relevant metadata
    with open(cdl_filename, "w", encoding='utf-8') as modified_cdl:
        modified_cdl.write(modified_cdl_content)

    return cdl_filename, nc_filename

def _cleanup_temporary_cdl(cdl_filename: Path) -> None:
    """Delete temporary CDL file after it's been used to create output netCDF

    :param cdl_filename: path to temporary CDL we just created
    :type cdl_filename: pathlib.Path
    :returns: {None}
    :rtype: {None}
    """

    logging.info('--- %-30s: %s', 'Cleaning up temporary CDL', cdl_filename)

    # delete temporary CDL, once its purpose is fulfilled
    os.remove(cdl_filename)

def _write_netcdf_from_cdl(
    cdl_filename: Path,
    nc_filename: Path,
    cleanup_cdl: bool = False) -> Dataset:
    """Write empty netCDF file using temporary CDL file we just created from template

    This function leverages ncgen to create an empty netCDF file which has structure, variables,
    and attributes as prescribed by temporary CDL file we just created from template. The actual
    ncgen command looks like: `ncgen -knc4 -o <path_to_netcdf> <path_to_cdl>`. This ncgen command
    is run via netCDF4.Dataset.fromcdl(). It also returns a netCDF4 Dataset object, which is used
    to populate the netCDF file with our data.

    Note that we write all netCDF output files to /data/output/... by default. To access this file
    elsewhere, one can mount a local directory into the container as /data/output or write the
    output file to s3 by providing a valid s3 URI to the -o command line argument.

    :param cdl_filename: path to temporary CDL we just created
    :type cdl_filename: pathlib.Path
    :param nc_filename: path to temporary CDL we just created
    :type nc_filename: pathlib.Path
    :param cleanup_cdl: whether or not we want to delete the temporary CDL after it's been used
    :type cleanup_cdl: bool
    :returns: Dataset object representation of netCDF file, used to populate file with our data
    :rtype: netcdf4.Dataset
    """

    logging.info('--- %-30s: %s', 'Writing netCDF from CDL', nc_filename)

    # create empty netCDF that has structure prescribed by temporary CDL
    output_dataset = Dataset.fromcdl(str(cdl_filename), str(nc_filename))

    # if we want to delete the temporary CDL, do so. Sometimes we don't for debugging purposes
    if cleanup_cdl:
        _cleanup_temporary_cdl(cdl_filename)

    return output_dataset

def prep_output_file(
    config_vars: dict[str, Any],
    interpolated_data: dict[str, xr.Dataset],
    cleanup_cdl: bool = False) -> tuple[Dataset, Path]:
    """Prepare netCDF file from CDL template and return netCDF object to populate

    This includes:
    - Writing temporary CDL file from template CDL file
    - Creating netCDF file from temporary CDL file (and returning Dataset object to populate)

    Note that we write all netCDF output files to /data/output/... by default. To access this file
    elsewhere, one can mount a local directory into the container as /data/output or write the
    output file to s3 by providing a valid s3 URI to the -o command line argument.

    :param config_vars: dict containing all attributes from ProcConfig class
    :type config_vars: dict
    :param interpolated_data: dict containing vertical coord-wise Xarray Datasets (interpolated)
    :type interpolated_data: dict
    :param cleanup_cdl: whether or not we want to delete the temporary CDL after it's been used
    :type cleanup_cdl: bool
    :returns: tuple containing Dataset object representation of netCDF file and filename
    :rtype: (netcdf4.Dataset, str)
    """

    # create temporary CDL file from template CDL file
    cdl_filename, nc_filename = _write_temporary_cdl(config_vars, interpolated_data)

    return _write_netcdf_from_cdl(cdl_filename, nc_filename, cleanup_cdl), nc_filename
