"""
Description
-----------

This workflow constitutes a Python port of Dave Turner's column extraction workflow, which
extracts "vertical profiles" of select variables from RAP/HRRR/GFS model output over sites
of interest. For each site, the nearest model point is identified and then bilinearly interpolated
in the horizontal and vertical (using geometric heights defined in configuration file). Data are
output to netCDF files.

Command Example
---------------

.. code-block:: console

    python3 src/main.py -f /data/input/<file name> -c ./conf
    -e <optional Minio s3 endpoint URL> -o <optional s3 URI>
    -m <optional model if file name doesn't specify> -l /data/log -g %Y%j%H%M%S -p ColumnExtractor

Command Line Arguments
----------------------

- -f, --fileToProc :pathlib.Path
    Path to file to process and extract columns from
    (required)
- -o, --outputDataDir :pathlib.Path
    Path to directory where processed output should land
    (required)
- -e, --endpointUrl :str
    Minio s3 endpoint URL used to fetch input files from internal /public
    (optional)
- -c, --configDir :pathlib.Path
    Path to configuration directory containing config files, profiling sites file,
    grib variable mapping, RAP grid, etc.
    (required)
- -l, --logDir :pathlib.Path
    Path to where log files should be written
    (optional; default is stdout)
- -g, --globPattern :str
    Pythonic format of job time string that is tacked on to the end of the log files
    (optional; default is ``%Y%j%H%M%S``)
- -p, --procName :str
    The process name
    (optional; default is ``ColumnExtractor``)

"""

from column_extractor.model_data import ModelData
from column_extractor.file_io import push_output_to_s3

def main() -> None:
    """Main method responsible for executing all phases of the workflow

    Workflow flows as follows:

    #. Initialize workflow configuration: parse and log command line arguments, set up logging,
       read configuration file, profiling sites, and grib variable mapping + extract for model.
    #. Preprocess data: extract all variables from file (as prescribed by grib variable mapping)
       and merge like vertical coordinates into single Xarray Dataset, doing all necessary
       massaging.
    #. Extract model data in NxN neighborhood surrounding closest model gridpoint to profiling site
    #. Post-process subsetted data: Recast winds as earth-relative, resample hybrid levels to
       uniform vertical levels, derive additional fields, perform unit conversions.
    #. Interpolate data to profiling sites.
    #. Write data at closest model gridpoints and interpolated to profiling sites to netCDF

    :param {None}: {None}
    :returns: {None}
    :rtype: {None}
    """

    # initialize instance of ModelData class and assign workflow configuration as attribute
    model_data = ModelData()

    # log command line arguments and model specific configuration parameters
    model_data.log_model_specific_configs()

    # pre-process input file
    preprocessed_data = model_data.preprocess()

    # for each profiling site of interest, identify closest model gridpoint
    # extract model data at each of these closest model gridpoints
    closest_gridpoint_data = model_data.extract_closest_gridpoint_data(preprocessed_data)

    # post-process model data (perform unit conversions, recast winds as earth-relative,
    # recast hybrid levels, derive additional fields, etc.)
    post_processed_data = model_data.post_process(closest_gridpoint_data)

    # interpolate model data to profiling sites of interest
    interpolated_data = model_data.interpolate_to_profiling_sites(post_processed_data)

    # write data at closest model gridpoint and interpolated to profiling site to netCDF
    local_filename = model_data.write_output_to_netcdf(post_processed_data, interpolated_data)

    # if an alternative output path is specified that references 's3', push output file to S3
    # (whether AWS or /public via Minio gateway)
    if model_data.input_output['output_path'] and 's3' in model_data.input_output['output_path']:
        push_output_to_s3(
            local_filename,
            model_data.input_output['output_path'],
            model_data.input_output['endpoint_url'])

    # log conclusion of job
    model_data.log_workflow_completion()

# run main routine
if __name__ == '__main__':
    main()
