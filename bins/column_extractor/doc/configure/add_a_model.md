# Adding a model

We often receive requests to add processing support for a new model

This is the most complex configuration change, as it requires:

- Updating all configuration files
- Wrangling with a new set of GRIB variables that may or may not
  play nice
- (potentially) extending the source code

Once added, follow the [build/deploy instructions](../build_and_deploy.md).

## How do I do it?

Say we want to add processing support for the RRFS-A experimental runs at NCEP.

First, we must choose a name for the model. This name is used in all configuration
files and is referenced throughout the application.

> [!CAUTION]
> Since the application attempts to determine which of the supported models
> an input file comes from, it's important that our choice matches the name of the model
> in their file names.

Let's choose `rrfs` to refer to the RRFS-A runs, since the file paths look like
`rrfs.tCCz.{natlev,prslev}.fHHH.grib2`.

### Add `rrfs` section to [`column_extractor.cfg`](../../conf/column_extractor.cfg)

These sections rarely differ (unless we need to provide the grid), so just copy and
paste another section and make sure `vertical_coord` is correct:

```python
[rrfs]
vertical_coord = hybrid
new_heights_m = 10,20,40,60,80,100,120,140,160,180,200,225,250,275,300,350,400,450,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,4200,4400,4600,4800,5000,5500,6000,6500,7000,8000,9000,10000,11000,12000,13000,14000,15000
neighbors = 1
```

### Add `rrfs` section to [`profiling_sites.json`](../../conf/profiling_sites.json)

As a first pass, this can just be a copy-and-paste of another model's sites (HRRR is
usually a good start). Request guidance from the requestor of this model for what
locations they'd like included.

> [!CAUTION]
> Ensure that all sites of interest are contained within the model domain.
> The application will drop sites during processing if found to be outside, or on the edge,
> of the model domain.

### Add `rrfs` section to [`variable_mapping.json`](../../conf/variable_mapping.json)

The specifics will depend on which variables are requested, which variables are available,
and how well they're handled by eccodes.

A combination of eccode's `grib_ls`/`grib_dump`, NCEP's `.idx` files, and
[ECMWF's GRIB parameter database](https://codes.ecmwf.int/grib/param-db/)
should get you pretty far. Use other model sections for pointers.
It's ultimately trial and error.

For extending supplementary GRIB tables as needed, see
[directions](../debug.md#fixing-unknown-grib-variables).

### Add `columns_output_rrfs_template.cdl`

Starting from another model's CDL file is a good start. The specifics depend on which
variables are requested and which variables are available. Some new variables
may need added, some unavailable variables may need removed. For variables available across multiple models, maintain the same variable naming convention.

> [!NOTE]
> It's important that each model CDL file is specific to that model. Update
> global attributes accordingly. Ensure that all variables listed are available for this
> model. Don't list variables that aren't available, as this leads to unecessary bloat.

### Extend source code (as needed)

The source code may need extended to various degrees to support a new model. Below are a few
examples, listed most common to least common:

- The list of supported models is hard-coded into `proc_config.py`
(ie `['hrrr', 'rap', 'gfs']`) in two places to ensure that the provided input file(s)
are associated with a model we can process. Add `rrfs` to these lists
- Various script/class/function docStrings list the models. Add `RRFS-A` here (not critical,
  but is nice)
- A rotation grid is needed to convert winds from grid-relative to earth-relative.
  This rotation grid is different for each model because it depends on the model's grid.
  A `rrfs` case is needed in `add_rot_grid()` of `model_utils.py` and a derivation function
  is needed in `calculations.py`. **We don't know what the specifics of this rotation
  calculation are at this time. We're trying to get that from the RRFS folks**
- Depending on the model, other pieces of the source code may need
  changed/extended/refactored. This is model-dependent.

## After the fact

We can compare output from this version of the application for models that we already
support (ie HRRR, RAP, GFS). However, we can't compare output from this version
of the application for RRFS, because the previous production version didn't support RRFS.

Now that the production version supports RRFS, we need to add RRFS testing to our output
testing suite.

Update the for loop in [`test_build.sh`](../../integration/test_build.sh) to include
`"rrfs"`:

```bash
for MODEL in "rap" "hrrr" "gfs" "rrfs"
```

Also add the requisite RRFS-A test data to our test dataset under
`/public/retro`.
