# Debugging the application

The ability to rapidly debug the application is helpful during
development. There are a couple ways to do this:

- Running the container image in interactive mode, allowing you to enter
  it via an interactive shell. Optionally also mounting your local `/src`
  and `/conf` directories into the container.
- Running our output testing suite to ensure your changes don't result
  in drift in the output, relative to our current production version.

## Running the container in "interactive" mode

Say you're building `v3.1.0` and want to debug it. You have a development `v3.1.0` container image, but continue to rapidly iterate on the code and need to test it.

Instead of debugging/running the code locally, you can use the container image as a testing environment by launching it in interactive mode:

```shell
docker run --rm -it --name column-extractor -v ./src:/app/src -v ./conf:/app/conf --entrypoint /bin/bash column_extractor:v3.1.0
```

This method overrides the Python command defined by `CMD` in the Dockerfile and opens an
interactive terminal session running a bash shell.

It also mounts your local `./src` and `./conf` directories into the
running container, overriding whatever code/configs that were built
with the image to use the current version of the code on your machine.
This prevents you from needing to rebuild the container image every time
you change the code during development. You can then run the Python application
manually, as described [here](run.md#running-the-container-in-interactive-mode).

If you mount your data directories via `-v ~/data:/data`, you can also use interactive mode to interrogate input GRIB files. This is valuable when troubleshooting configuration / adding new models. eccodes has a variety of [CLI tools](https://confluence.ecmwf.int/display/ECC/GRIB+tools)
that can be used to interrogate GRIB files, namely `grib_ls` and
`grib_dump`.

**By debugging within a container instance instead of locally, we ensure that our debugging is
self-contained and not impacted by differences in our development machines.**

## Test a container build

Currently, the only tests we perform are on the output produced by different versions of the code. This is done to ensure that new versions of the code will produce identical output as the current production version (when applicable).

It's important to compare output for all models (HRRR, RAP, GFS) and for
analysis (ie FH00) and forecast (ie FH01+).
[test_build.sh](../integration/test_build.sh) will run the production
container image (`-t`) and development container image (`-d`) for each model
and two forecast hours (00, 01+) and compare the output using
[test_output.py](../src/tests/test_output.py).

[test_output.py](../src/tests/test_output.py) tests the following:

- Dimensions and their values match
- Minimum metadata attributes are present (though some of their values may differ)
- All variables exist and their dimensions/attributes/values match (we
  perform exact matches on variable values using Numpy's
  [`assert_array_equal()`](https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_array_equal.html))

To test output drift between production `v3.0.0` and development `v3.1.0`, run:

```shell
bash integration/test_build.sh -t v3.0.0 -d v3.1.0 &> development/v3.0.0_v3.1.0.txt
```

This will grab test input files for RAP, HRRR, and GFS that we have stored under
`/public/retro` and run the production/development containers on them as
discussed above.

> [!NOTE]
> If you already have the necessary test data under `./development/data/input`,
> you can use the `-l` flag to skip retrieving the files from `/public/retro`
> and use your local files.

> [!NOTE]
> The production/development container images must both be readily available
> locally in order to run them and then test their output.

For each model / time step, successful output from
[test_output.py](../src/tests/test_output.py) will look like:

```
test_output_metadata (tests.test_output.TestOutput)
Is minimum file metadata present? ... ok
test_output_structure (tests.test_output.TestOutput)
Is the output structure/dimensions unchanged? ... ok
test_variables (tests.test_output.TestOutput)
Are variables present and do their dimensions/attributes/content match? ... ok

----------------------------------------------------------------------
Ran 3 tests in 9.116s

OK
```

## Fixing unknown GRIB variables

There are instances where the GRIB encoding for certain variables isn't recognized
by eccodes. This may be due to the GRIB table version or (oftentimes) differences between
NCEP GRIB conventions and ECMWF GRIB conventions (for which eccodes was built).

eccodes assigns names/units/etc. of `unknown` to variables it doesn't recognize. Sometimes
these unknown variables are ones that we want. In order to make them explicit, we need
to add these variables to the supplementary GRIB definitions, located under
[`conf/custom_definitions/grib2/localConcepts`](../conf/custom_definitions/grib2/localConcepts).

> [!NOTE]
> eccodes first checks it's built-in GRIB tables, then checks the supplementary GRIB
> definitions under [`conf/custom_definitions/grib2/localConcepts`](../conf/custom_definitions/grib2/localConcepts).
> NCEP is referred to as `kwbc`, so NCEP custom definitions sit under
> the [`kwbc subdirectory`](../conf/custom_definitions/grib2/localConcepts/kwbc).
> For GSL experimental runs, the center may be different, so a new subdirectory may
> be needed.

`conf/custom_definitions/grib2/localConcepts/{center}` contains several files that define
a variable's name, description, units, etc. Each of these attributes is tied together
by the variable's `discipline`, `parameterCategory`, and `parameterNumber`. EG:

```
# Planetary boundary layer height
'hpbl' = {
         discipline = 0 ;
         parameterCategory = 3 ;
         parameterNumber = 196 ;
}
```

Once you know the `discipline`, `parameterCategory`, and `parameterNumber` used to
encode the unknown GRIB variable you're interested in, you can add the corresponding
attributes to each of the files under `conf/custom_definitions/grib2/localConcepts/kwbc`. Of particular importance are the
variable's `shortName` and `units`, because these are used within
the application.

Occasionally, you may find that a vertical coordinate type is unknown. In this case, a
variable might be decoded correctly but `typeOfLevel` will be `unknown`. Identify the corresponding `typeOfFirstFixedSurface`
and add it to [`typeOfLevelConcept.def`](../conf/custom_definitions/grib2/localConcepts/kwbc/typeOfLevelConcept.def).

`discipline`, `parameterCategory`, `parameterNumber`, and `typeOfLevelConcept` can all
be found via eccode's `grib_dump` and by discussing with the model output file producers
(when possible).