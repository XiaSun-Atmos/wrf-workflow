"""
Description
-----------

Class TestOutput runs a variety of unit tests to compare a "truth" output file
against an output file from code that's in development. A variety of components are
verified between the truth and development output, including:

- Same file dimensions (by name)
- Each dimension has the same size and values (except for profiling sites)
- Minimum metadata present as global attributes
- Each variable from truth output is present
- Variable dimensions are unchanged
- Variable attributes are unchanged
- Variable is within 4 decimal places of truth

"""

import unittest
import logging
from pathlib import Path
import xarray as xr
import numpy as np

class TestOutput(unittest.TestCase):
    """Runs a variety of unit tests to compare output file contents

    This class compares the dimensionality, metadata, and variables of output files from the
    current production code (ie "truth") and the in-development code to verify that output
    hasn't drifted for unknown reasons.

    """

    @classmethod
    def setUpClass(cls):
        """Read in "truth" output and new output to compare against. Runs at TestOutput init"""

        # define truth and test output files
        truth_output_file = Path('/data/testing/truth_output.nc')
        test_output_file = Path('/data/testing/test_output.nc')

        # check if truth/test output files exist
        for file in [truth_output_file, test_output_file]:
            if not file.is_file():
                raise FileNotFoundError(f'{file} not found, exiting...')

        # read in test/truth output
        truth_output = xr.open_dataset(truth_output_file, engine='netcdf4')
        test_output = xr.open_dataset(test_output_file, engine='netcdf4')

        # assign coordinates profiling site / height / depth / forecast hour coordinates
        truth_output = truth_output.assign_coords(
            {'sites':truth_output['dsite'],
            'height_bin':truth_output['height'],
            'depth_bin':truth_output['depth'],
            'forecast_hour':truth_output['forecast']})
        # drop the variables we assigned as coordinates
        truth_output = truth_output.drop_vars(['dsite','height','depth','forecast'])

        # assign coordinates profiling site / height / depth / forecast hour coordinates
        test_output = test_output.assign_coords(
            {'sites':test_output['dsite'],
            'height_bin':test_output['height'],
            'depth_bin':test_output['depth'],
            'forecast_hour':test_output['forecast']})
        # drop the variables we assigned as coordinates
        test_output = test_output.drop_vars(['dsite','height','depth','forecast'])

        # assign as class attributes
        cls.truth_output = truth_output
        cls.test_output = test_output

    def test_output_structure(self):
        """Is the output structure/dimensions unchanged?"""

        # first test to make sure dimensions are the same
        self.assertEqual(
            TestOutput.truth_output.sizes.keys(),
            TestOutput.test_output.sizes.keys(),
            "Dimensions are different!")

        # then test to make sure dimensionality is the same
        for truth_dim in TestOutput.truth_output.dims:
            with self.subTest(truth_dim):
                if truth_dim == 'sites':
                    # don't test profiling sites, just warn if different
                    # (since these change so often)
                    if (TestOutput.truth_output.sizes[truth_dim] !=
                        TestOutput.test_output.sizes[truth_dim]):
                        logging.warning('# of profiling sites has changed!')
                else:
                    np.testing.assert_array_equal(
                        TestOutput.truth_output[truth_dim].data,
                        TestOutput.test_output[truth_dim].data,
                        f'{truth_dim} dimension has changed!')

    def test_output_metadata(self):
        """Is minimum file metadata present?"""

        # loop through each attribute in truth output (considered the "minimum" metadata)
        # and make sure it's present in new output
        for truth_attr in TestOutput.truth_output.attrs:
            with self.subTest(truth_attr):
                self.assertIn(
                    truth_attr,
                    TestOutput.test_output.attrs,
                    f'{truth_attr} metadata is missing!')

    def test_variables(self):
        """Are variables present and do their dimensions/attributes/content match?"""

        # loop through each variable in truth output and check for:
        # 1) variable exists
        # 2) variable dimensions match
        # 3) variable attributes match
        # 4) variable content matches

        # test model initialization times first, then drop
        with self.subTest('model_initialization'):
            np.testing.assert_array_equal(
                TestOutput.truth_output['model_initialization'].data,
                TestOutput.test_output['model_initialization'].data,
                'Model initialization times differ!')

        # drop model initialization times
        truth_output = TestOutput.truth_output.drop_vars(['model_initialization'])
        test_output = TestOutput.test_output.drop_vars(['model_initialization'])

        # first, find profiling sites shared by truth and test output
        # only compare variable values between shared sites
        shared_sites = np.intersect1d(
            truth_output['sites'].values,
            test_output['sites'].values)

        # test all site-wise variables
        for truth_var_name,truth_var_data in truth_output.data_vars.items():
            with self.subTest(truth_var_name):
                # check that variable exists in test output
                self.assertIn(
                    truth_var_name,
                    test_output.data_vars,
                    f'{truth_var_name} is missing!')
                # check that variable dimensions match those in test output
                self.assertEqual(
                    truth_var_data.dims,
                    test_output[truth_var_name].dims,
                    f'{truth_var_name} dimensions are different!')
                # check that variable attributes match those in test output
                self.assertEqual(
                    truth_var_data.attrs,
                    test_output[truth_var_name].attrs,
                    f'{truth_var_name} attributes are different!')
                # loop through each shared site and compare
                for site in shared_sites:
                    with self.subTest(site):
                        np.testing.assert_array_equal(
                            truth_var_data.sel(sites=site).data,
                            test_output[truth_var_name].sel(sites=site).data,
                            f'Values in {truth_var_name} have changed at {site}!')

if __name__ == '__main__':
    unittest.main()
