"""
Description
-----------

Class ProcConfig handles workflow initialization tasks (argument parsing, log handling,
reading configuration files, retrieving input file if needed, etc.). It subsets the
configuration files by the model that the input file is associated with (ie RAP, HRRR, GFS, etc).
These model-specific parameters are then set as attributes of the class.

"""

import argparse
import sys
import logging
import time
import configparser
from typing import Optional
from pathlib import Path
from datetime import datetime
from column_extractor import file_io

class ProcConfig:
    """Handles workflow initialization tasks and model-specific configuration tasks

    Attributes:
        start_time (str): The start time of the workflow.
        input_output (dict): File to process (local/remote path), endpoint URL, output path,
        path to template CDL for model
        model (str): The supported model associated with the file.
        model_config (dict): Model-specific configuration parameters.
        model_profiling_sites (dict): Model-specific profiling sites.
        model_var_map (dict): Model-specific GRIB variable mapping.
    Methods:
        log_model_specific_configs(): Logs the model-specific configuration parameters, profiling
            sites, and variables to extract.
        log_workflow_completion(): Upon completion of workflow, log this + time elapsed

    Usage:
        # Create an instance of ProcConfig
        config = ProcConfig()

        # Access the attributes and methods of the instance
        print(config.<attribute>)
        config.<method>()
    """

    def __init__(self) -> None:
        """Initialize the ProcConfig instance

        Performs startup tasks such as parsing command line arguments, retrieving configuration
        contents, reading the profiling sites file, and obtaining the GRIB variable mapping.

        :param {None}: {None}
        :returns: {None}
        :rtype: {None}
        """

        # perform startup tasks (parse command line arguments, retrieve configuration
        # contents, read profiling sites file, grib variable mapping file)

        # document start time of workflow
        self.start_time = datetime.utcnow()

        # process arguments from the command line
        args = self._get_input_args()

        # set up logger
        self._set_up_logging_logistics(args)

        # log workflow start
        logging.info('Starting workflow at %s UTC...\n',
            self.start_time.strftime("%Y-%m-%d %H:%M:%S"))

        # log command line arguments
        self._log_input_params(args)

        # process file input/output parameters; fetch from AWS/public via S3 protocol if needed
        self._get_set_input_output(args)

        # if file(s) to process doesn't exist, raise FileNotFoundError
        for file in self.input_output['local_paths']:
            file_io.file_exists(file)

        # read in main configuration file ('column_extractor.cfg')
        config_contents = file_io.read_config(args.config_dir / 'column_extractor.cfg')

        # read in profiling sites file ('profiling_sites.json')
        all_profiling_sites = file_io.read_json(args.config_dir / 'profiling_sites.json')

        # read in GRIB variable mapping file ('variable_mapping.json')
        grib_var_map = file_io.read_json(args.config_dir / 'variable_mapping.json')

        # set ProcConfig model-specific attributes depending on inputs
        self._set_model(args.user_defined_model, args.file_locations)
        self._set_model_config(config_contents)
        self._set_model_profiling_sites(all_profiling_sites)
        self._set_model_var_map(grib_var_map)
        self._get_set_model_output_template(args.config_dir)

        # verify that we have minimum necessary input files based on file types
        # referenced in model variable map
        self._verify_minimum_files()

    @staticmethod
    def _get_input_args() -> argparse.Namespace:
        """Parse arguments passed to the program via the command line

        :param {None}: {None}
        :returns: parsed command line arguments
        :rtype: argparse.Namespace
        """

        # establish argument parser object
        parser = argparse.ArgumentParser(
            description='Extract vertical profiles (AKA columns) over locations of note from'
            ' GRIB format model grids and output as netCDF')

        # add input arguments to the argument parser object
        parser.add_argument(
            '-m',
            '--userDefinedModel',
            dest='user_defined_model',
            type=str,
            choices=['rap', 'hrrr', 'gfs', 'wfip1km', 'wfip3km'],
            help='model that file is associated (ie RAP, HRRR, GFS, etc.)')
        parser.add_argument(
            '-f',
            '--fileLocations',
            dest='file_locations',
            type=str,
            nargs='+',
            required=True,
            help='original location of file(s) to be processed')
        parser.add_argument(
            '-e',
            '--endpointUrl',
            dest='endpoint_url',
            type=lambda value: None if value in ('None', 'none') else value,
            help='Optional s3 endpoint URL for /public access via Minio')
        parser.add_argument(
            '-c',
            '--configDir',
            dest='config_dir',
            type=Path,
            required=True,
            help='parent directory for all configuration-type files')
        parser.add_argument(
            '-l',
            '--logDir',
            dest='log_dir',
            type=lambda value: None if value in ('None', 'none') else Path(value),
            help='parent directory path where log files should be written')
        parser.add_argument(
            '-o',
            '--outputPath',
            dest='output_path',
            type=lambda value: None if value in ('None', 'none') else value,
            help='destination to write output file to (only s3 locations currently supported)')
        parser.add_argument(
            '-g',
            '--globPattern',
            dest='glob_pattern',
            default='%Y%j%H%M%S',
            help='pythonic format of job time string that is tacked on to the end of state file')
        parser.add_argument(
            '-p',
            '--procName',
            dest='proc_name',
            default='ColumnExtractor',
            help='the process name (default is ColumnExtractor)')

        # extract arguments from the argument parser object and return argparse.Namespace object
        return parser.parse_args()

    def _set_up_logging_logistics(self, args: argparse.Namespace) -> None:
        """Set up logging: log file syntax, timestamp convention, etc.

        If log directory specified, log files will be written there. If no log directory
        specified, logs will be written to stdout

        :param args: command line arguments
        :type args: argparse.Namespace
        :returns: {None}
        :rtype: {None}
        """

        # if log directory specified, construct logging handler here
        # log files are job-wise, with job time used to differentiate
        if args.log_dir:
            if not args.log_dir.is_dir():
                args.log_dir.mkdir(parents=True)
            log_file = args.log_dir / (f'{args.proc_name}.log'
                f'.{self.start_time.strftime(args.glob_pattern)}')
            handler = logging.FileHandler(filename = log_file)
            # declare basic config
            logging.basicConfig(
                handlers=[handler],
                format='%(asctime)s %(levelname)s: %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S',
                level=logging.INFO)
        # if no log directory specified, then configure to print to stdout
        else:
            logging.basicConfig(
                stream=sys.stdout,
                format='%(asctime)s %(levelname)s: %(message)s',
                datefmt='%Y-%m-%d %H:%M:%S',
                level=logging.INFO)

        logging.Formatter.converter = time.gmtime

    @staticmethod
    def _log_input_params(args: argparse.Namespace) -> None:
        """Log input parameters as prescribed by ``args``

        :param args: command line arguments
        :type args: argparse.Namespace
        :returns: {None}
        :rtype: {None}
        """

        # log command line arguments
        logging.info('Command line arguments are:')
        for arg, value in vars(args).items():
            logging.info('--- %-25s: %s', arg, value)
        logging.info('Moving on...\n')

    def _get_set_input_output(self, args: argparse.Namespace) -> None:
        """Get/set file input/output parameters

        If file location(s) includes reference to s3, attempt to download file(s) via s3 protocol
        to /data/input. Set input_output dict with remote path(s) (if applicable), local path(s),
        endpoint url (if applicable), and output path(s).

        :param args: command line arguments
        :type args: argparse.Namespace
        :returns: {None}
        :rtype: {None}
        """

        # instantiate empty lists for local/remote path(s)
        local_paths = []
        remote_paths = []

        for file in args.file_locations:
            # check to see if files exist somewhere in s3 (whether AWS s3 bucket
            # or s3 connector to /public via Minio)
            if 's3' in file:
                # if so, check to see if it's a correctly formatted S3 path
                if parsed_uri := file_io.validate_s3_uri(file, args.endpoint_url):
                    # if valid s3 URI, store remote path
                    remote_paths.append(file)
                    # if valid s3 URI, attempt to download file
                    local_paths.append(file_io.download_from_s3(parsed_uri, args.endpoint_url))
            # if files contains no reference to an s3 path, assume that the file is
            # locally available and proceed as normal
            else:
                logging.info('File (should be) available locally')
                logging.info('--- %-20s: %s\n', 'Path is', file)
                # remote path doesn't exist, set to None
                remote_paths.append(None)
                # store local path
                local_paths.append(Path(file))

        # input_output attribute contains local path, remote path, endpoint url (if applicable),
        # and output path
        self.input_output = {
            'remote_paths': remote_paths,
            'endpoint_url': args.endpoint_url,
            'local_paths': local_paths,
            'output_path': args.output_path}

    def _set_model(self, user_defined_model: Optional[str], file_locations: list) -> None:
        """Set the model attribute based on the input file(s) OR command line argument

        If model is defined via the command line, use that. If not, test input file(s) against our
        currently supported models (via string comprehension). If model not specified in input
        filename(s), then must use command line argument. If we don't support the prescribed model
        or if files come from multiple models, the program raises a RuntimeError.

        :param user_defined_model: model that file comes from, defined via command line
        :type user_defined_model: str
        :param file_locations: where the model file(s) are located
        :type file_locations: list
        :returns: {None}
        :rtype: {None}
        """

        # if user already specified model via command line argument, use it
        if user_defined_model:
            self.model = user_defined_model
        else:
            # check if input filename(s) contain reference to any of our currently supported models
            models = []
            for file in file_locations:
                # loop through supported models, check if referenced in input filename
                # if it is, store. if no match, store NOT_SUPPORTED
                models.append(
                    next((model for model in ['hrrr', 'rap', 'gfs', 'wfip1km', 'wfip3km'] if model in file.lower()),
                    'NOT_SUPPORTED'))
            # if NOT_SUPPORTED for any input file(s), raise RuntimeError
            if 'NOT_SUPPORTED' in models:
                raise RuntimeError('One or more unknown models. Exiting...\n')
            # if more than one unique model, raise RuntimeError
            if len(set(models)) > 1:
                raise RuntimeError('Multiple models found. Exiting...\n')
            # otherwise, set model attribute using first entry (since all entries are the same)
            self.model = models[0]

    def _set_model_config(self, config_contents: configparser.ConfigParser) -> None:
        """Set the model_config attribute based on the provided configuration contents

        Retrieves model-specific configuration contents from the config_contents dictionary.

        :param config_contents: contents of configuration file (all models)
        :type config_contents: configparser.ConfigParser
        :returns: {None}
        :rtype: {None}
        """

        # retrieve model-specific config contents
        self.model_config = config_contents[self.model]

    def _set_model_profiling_sites(self, all_profiling_sites: dict) -> None:
        """Set the model_profiling_sites attribute based on the provided profiling sites

        Retrieves model-specific profiling sites from the all_profiling_sites dictionary.

        :param all_profiling_sites: contents of profiling sites file (all models)
        :type all_profiling_sites: dict
        :returns: {None}
        :rtype: {None}
        """

        # retrieve model-specific profiling sites
        self.model_profiling_sites =  all_profiling_sites[self.model]

    def _set_model_var_map(self, var_map: dict) -> None:
        """Set the model_var_map attribute based on the provided GRIB variable mapping

        Retrieves model-specific GRIB variable mapping from the var_map dictionary.

        :param var_map: contents of grib variable mapping file (all models)
        :type var_map: dict
        :returns: {None}
        :rtype: {None}
        """

        # retrieve model-specific grib variable mapping
        self.model_var_map =  var_map[self.model]

    def _get_set_model_output_template(self, config_directory: Path) -> None:
        """Get/set the model_output_template attribute based on the model

        Checks for existence of model CDL file from which output file is created.

        :param config_directory: path to configuration directory
        :type config_directory: pathlib.Path
        :returns: {None}
        :rtype: {None}
        """

        # construct model output template path
        model_output_template_path = config_directory / f'columns_output_{self.model}_template.cdl'

        # check existence
        file_io.file_exists(model_output_template_path)

        # assign to attribute
        self.input_output['model_output_template_path'] = model_output_template_path

    def _verify_minimum_files(self) -> None:
        """Verify that we have the minimum input files

        Each model variable mapping in variable_mapping.json lists a file type from which to
        retrieve each variable from. Here, we ensure that we have the minimum input files needed
        based on the referenced file types.

        :param {None}: {None}
        :returns: {None}
        :rtype: {None}
        """

        # compile all needed file types from model variable mapping pre-process section
        file_types = []
        for coord_vars in self.model_var_map['pre-process'].values():
            for coord_var in coord_vars:
                file_types += coord_var['file_id']

        # verify that we have an input file corresponding to each file type that we need
        have_all_file_types = True
        for file_type in set(file_types):
            if not any(file_type in str(file) for file in self.input_output['local_paths']):
                have_all_file_types = False

        # if we have all necessary files, log and proceed
        if have_all_file_types:
            logging.info('Have all specified input file types. Proceeding...\n')
        # if we don't, raise RuntimeError
        else:
            raise RuntimeError('One or more necessary file types not present. Exiting...\n')

    def log_model_specific_configs(self) -> None:
        """Log model-specific configuration parameters, profiling sites, and variables to extract

        :param {None}: {None}
        :returns: {None}
        :rtype: {None}
        """

        #log model
        logging.info('Model is: %s\n', self.model.upper())

        # log config parameters
        logging.info('%s configuration parameters are:', self.model.upper())
        for key in self.model_config:
            logging.info('--- %-20s: %s', key, self.model_config[key])
        logging.info('Moving on...\n')

        # log profiling sites
        logging.info('%s profiling sites are:', self.model.upper())
        for site_name, site_coords in self.model_profiling_sites.items():
            logging.info('--- %-70s: %-10s, %+10s', site_name, site_coords[0], site_coords[1])
        logging.info('Moving on...\n')

        # log grib variables
        logging.info('%s variables to extract are:', self.model.upper())
        for vertical_coord, vars_ in self.model_var_map['pre-process'].items():
            for var_info in vars_:
                filter_key_info = [
                    f"{key+':':<15s}{value:<20s}" for key,value in var_info['filter_keys'].items()]
                logging.info('--- %s%s',
                    f"{'typeOfLevel:':<20s}{vertical_coord:<25s}",
                    ''.join(filter_key_info))
        logging.info('Moving on...\n')

    def log_workflow_completion(self) -> None:
        """Upon completion of workflow, log this + time elapsed

        :param {None}: {None}
        :returns: {None}
        :rtype: {None}
        """

        # log conclusion of job
        logging.info(
            '%-25s: %s UTC',
            'Workflow Completed At', datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
        logging.info('%-25s: %s', 'Time Elapsed', datetime.utcnow() - self.start_time)
        logging.info('--------------------------------------------------------')
        logging.info('--------------------------------------------------------\n')
