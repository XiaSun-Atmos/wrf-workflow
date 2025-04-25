# Adding a variable

We may receive a request to add a certain variable to the output
for a model(s).

This isn't a request we often get, but is worth considering.

The ease of this change varies depending on the variable, model,
and whether or not GRIB tables need to be modified.

Once added, follow the [build/deploy instructions](../build_and_deploy.md).

## How do I do it?

To add a variable, you'll need to verify:

- That the model(s) output this variable in the first place
- Whether the requestor wants the instantaneous, accumulated, or cycle-average
  version of the variable (if applicable)
- Which model output file type(s) this variable is in (ie `wrfnat`, `wrfprs`, `pgrb2`, `pgrb2b`)
- The eccodes `typeOfLevel`, `shortName`, `stepType`, and `stepRange` (if applicable)
  for this variable
- An acceptable variable name for the output files

Say we want to output Convective Available Potential Energy (CAPE) for HRRR. We see that CAPE
is available in HRRR's `wrfnat` files, by interrogating an `.idx` file. We know it's an
instantaneous field because it's only valid for a single forecast hour.

Now we need to know how eccodes decodes it. It's easiest to use eccode's `grib_ls` CLI tool
to determine the correct `typeOfLevel`, `shortName`, `stepType`, and `stepRange` (if applicable).
This can be done from within the container in
[interactive mode](../run.md#running-the-container-in-interactive-mode).

In doing so, you'll arrive at:

```
(base) column_extractor@c7818fc09dfc:~$ grib_ls /data/columns/input/hrrr.t00z.wrfnatf13.grib2
/data/columns/input/hrrr.t00z.wrfnatf13.grib2
edition    centre    date    dataType    gridType    stepRange    typeOfLevel    level    shortName    packingType  
...
2    kwbc    20240508    fc    lambert    13    surface    0    cape    grid_complex_spatial_differencing 
```

Now we know the following:

- `shortName` is `cape`
- `typeOfLevel` is `surface`
- `stepType` is `instant`

> [!NOTE]
> You can use [ECMWF's GRIB parameter database](https://codes.ecmwf.int/grib/param-db/)
> to get more information on which variables correspond to which `shortName`. Note that
> not all of NCEP's GRIB codes will be included in this.

### Extend [`variable_mapping.json`](../../conf/variable_mapping.json)

We'll need to update the HRRR `pre-process` and `post-process` sections.

In the `surface` section of `pre-process`, you'll add the following:

```json
{
    "filter_keys": {
        "shortName": "cape",
        "stepType": "instant"
    },
    "file_id": ["wrfnat"]
}
```

In the `surface` section of `post-process`, you'll add the following:

```json
{
    "var_name": "cape",
    "netcdf_mapping": ["cape"]
}
```

> [!NOTE]
> The order of each vertical coordinate type's variables doesn't matter.

### Extend the relevant CDL file

We'll need to update
[`columns_output_hrrr_template.cdl`](../../conf/columns_output_hrrr_template.cdl) as well.

In the `variables` section, add:

```
float cape0(forecast_hour, sites, height_bin) ;
        cape0:long_name = "Convective available potential energy at the closest model grid location" ;
        cape0:units = "J/kg" ;
        cape0:_FillValue = 9.96921e+36f ;
float cape1(forecast_hour, sites, height_bin) ;
        cape1:long_name = "Convective available potential energy calculated by bilinear interpolation via Scipy's griddata()" ;
        cape1:units = "J/kg" ;
        cape1:_FillValue = 9.96921e+36f ;
```

> [!NOTE]
> The order of variables in the CDLs doesn't matter.

## I can't find the variable I'm after

If you've verified that the variable you want to add is output by the model(s)
in question but cannot find it in the `grib_ls` output, it's likely that
it's being decoded as `unknown`.

This means you'll need to add it to the supplementary GRIB definitions. See
[directions](../debug.md#fixing-unknown-grib-variables).
