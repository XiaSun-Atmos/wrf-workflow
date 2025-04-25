# Adding a site of interest

We often receive requests to extract data at additional sites of interest.
This is the least complicated configuration change to make.

## How do I do it?

The sites of interest at which we extract data are defined in
[`profiling_sites.json`](../../conf/profiling_sites.json).

To add a site, you'll need:

- A name for the site
- A latitude (in degrees N)
- A longitude (in degrees E)
- The model(s) it applies to

Once you have this information, add it to the relevant model
section(s) in `profiling_sites.json`.

Say we want to add Fort Collins, CO (40.5853, -105.0844) to HRRR.
To the `hrrr` section of `profiling_sites.json`, you'd add:

```json
"Fort Collins, CO": [
    40.5853,
    -105.0844
]
```

> [!NOTE]
> There are no particular conventions for site names and the sites
> aren't in any particular order. However, try to follow the
> naming conventions of like sites and try to group like sites
> together logically.

Once added, follow the [build/deploy instructions](../build_and_deploy.md).

## Automating via GitHub Actions

We've prototyped a GitHub Action that'll automatically add the requested
site of interest and open a PR.

This Action is triggered by the opening of a GitHub Issue via an
[issue template](../../.github/ISSUE_TEMPLATE/new_profiling_site.yml).
This template includes boxes for site name, latitude, longitude and a dropdown
to select models from. An issue can be opened via this template by navigating
to [New Issue](https://github.com/NOAA-GSL/column-extractor/issues/new/choose)
and clicking `Get Started` on `New Profiling Site`.

When an issue is opened, the Action parses the issue fields, checks out
a new branch, inserts the new site information into `profiling_sites.json`
via jq, commits it, and opens a PR to main.

The thought was that Dave and/or DSG could open an issue and have the
site added automatically, saving us some time. This could be further extended
with a CI/CD component.

> [!CAUTION]
> This has been demonstrated to work, but hasn't been rigorously tested.
> Of note, some sites have latitude/longitude ints represented as floats
> (ie 26.0) in `profiling_sites.json`. When adding a site, jq changes all
> of these floats to ints (ie 26.0 -> 26). This results in unintended changes
to `profiling_sites.json`. However, I'm not sure if this actually results
in meaningful differences within the application itself.