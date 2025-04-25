# Building & deploying the application

Say you've updated the source code, added a new model, tweaked the model configurations, etc.
These changes will become `v3.1.0`. Now what?

## Prerequisites to deploying

Before proceeding, make sure you've:
- Committed your changes on a separate branch and pushed them to the repository.
- Built a development container image (you can preemptively tag this image
  as `v3.1.0`). See [build instructions](#building-the-container-image-locally).
- Verified that your changes don't impact the output by running our test suite,
  using the current production version as the "truth". See
  [testing instructions](debug.md#test-a-container-build) for more info.
- Opened a PR, and had this PR approved/merged into the main branch by one of your peers.
- Tagged this new version of the code with an annotated tag using semantic versioning
  (ie `v3.1.0`). Pushed this tag to the repository.

Once the above are done, we can start the deployment process.

## Building the container image locally

Say you want to build `v3.1.0`:
```shell
docker build --no-cache --pull --build-arg="ARG_CODE_VERSION=v3.1.0" -t column_extractor:v3.1.0 .
```

> [!NOTE]
> It's important that `ARG_CODE_VERSION` matches the tag applied via `-t column_extractor:<tag>` and
> that both of these match the corresponding git tag, if applicable. `ARG_CODE_VERSION` is assigned
> as an environment variable at build time and is added to each output file as metadata at runtime.

## Pushing the container image to repositories

We store our container images in two repositories: GSL's development [Harbor registry](https://harbor-dev.gsd.esrl.noaa.gov/harbor/projects)
and AWS Elastic Container Registry (ECR).

The images stored in Harbor are mainly there for posterity's sake and are easily accessible from within
GSL's network.

The images stored in AWS ECR are what we use for the real-time processing.

### Staging a container image

It's important that we don't push a soon-to-be production image to AWS when a model cycle is still
being processed. That'd result in a cycle having output produced by two different versions of the image.

To avoid this, our real-time processing in AWS uses the container image tagged `PRODUCTION`. That way,
we can push `v3.1.0`, without it being used until we update the tags to also include `PRODUCTION`.

To push a local `v3.1.0` container image to Harbor (dev) and AWS ECR, navigate to the root
of the repository and run:

```shell
bash ./integration/build_tag_push.sh -v v3.1.0 -p gsl-prod
```

> [!NOTE]
> Make sure you've authenticated to AWS via `aws sso login --profile gsl-prod`

### "Graduating" a container image to `PRODUCTION`

Once you've verified that all forecast hours from current cycles have been processed, it's
time to add the `PRODUCTION` tag to the `v3.1.0` container image. This will designate it as
the current production image.

To verify that no tasks are in progress, navigate to the `dsg-data-processing` ECS cluster in
GSL's prod AWS environment, and verify that there are no active tasks running. Also verify that all output files have been published by checking bucket `gsl-its-dsg-columns` in prod. Cycles usually finish being processed in the latter half of the hour (unless it's an extended HRRR/RAP cycle).

Once verified, run:

```shell
bash ./integration/build_tag_push.sh -v v3.1.0 -t PRODUCTION -p gsl-prod
```

This will add the `PRODUCTION` tag to the `v3.1.0` container image that already exists in Harbor dev
and AWS ECR. The next set of ECS tasks will now use `v3.1.0` as the production container image.

## Creating a new release on GitHub

Now that we're running a new version of the code in production, you'll
[create a new GitHub release](https://github.com/NOAA-GSL/column-extractor/releases/new)
to track changes.

Use tag `v3.1.0` and generate the automatic release notes.

Then publish release.
