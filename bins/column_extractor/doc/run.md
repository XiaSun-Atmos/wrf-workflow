# Running the application

The container image can run anywhere that Docker is supported. This includes your local
work station, the cloud, VMs, etc. 

The application can be run in two ways:
- By running the container image as an executable and assigning environment variables to the
  container at runtime.
- By running the container image in interactive mode, opening a terminal session within, and
  running the application via a standard `python3 ...` call with command line arguments.

The map between Docker environment variables and Python command line arguments for
relevant input is outlined in the run command at the bottom of the
[Dockerfile](../Dockerfile):

```docker
CMD /opt/conda/bin/python3 \
    src/main.py \
    -c \
    ./conf \
    -f \
    $FILE_LOC \
    -e \
    $ENDPOINT_URL \
    -o \
    $OUTPUT_PATH \
    -l \
    $LOG_DIR
```

> [!NOTE]
> You must have access to the container image in order to run it. You can build a particular version
> of the application by checking out the relevant version tag from GitHub and following the
> [build instructions](build_and_deploy.md#building-the-container-image-locally). You can also run the image for any of the 5
> most recent versions by first pulling it from Harbor (dev) or AWS ECR.

## Running the container as executable

This method assigns environment variables to the container at runtime, which are then passed to the Python command that runs when the container runs. This tells the application what files to process, where to put logs, etc.

> [!NOTE]
> This is how the container is executed in AWS as part of the real-time processing pipeline. We
> pass input file paths, output path, etc. as environment variable overrides when we submit
> an ECS task.

### Fetch remote file via S3 and process locally

The application supports s3 protocol for fetching files from publically accessible s3 buckets (ie NODD) and GSL internal `/public` (via Minio). When fetching from internal `/public`, the endpoint
URL `https://minio1.gsd.esrl.noaa.gov:9001` must be passed via the `ENDPOINT_URL` environment variable.

Here's how you'd run a local build of version `v3.0.0` to process a specific forecast hour of HRRR
from the NODD:

```shell
docker run --rm --name column-extractor -v ~/output:/data/output -e FILE_LOC="s3://noaa-hrrr-bdp-pds/hrrr.20240701/conus/hrrr.t00z.wrfnatf00.grib2 s3://noaa-hrrr-bdp-pds/hrrr.20240701/conus/hrrr.t00z.wrfprsf00.grib2" -e ENDPOINT_URL=None -e LOG_DIR=/data/output -e OUTPUT_PATH=None column_extractor:v3.0.0
```

> [!NOTE]
> `-v` represents a mount. `-v ~/data/output:/data/output` mounts your local `~/data/output` to the 
> container's `/data/output` (where the container writes its output file). This is only needed
> if you wish to interrogate the output file locally after the container exits.
> If not included, the output file will disappear when the container exits.

To process files from internal `/public`, you'd pass `-e ENDPOINT_URL=https://minio1.gsd.esrl.noaa.gov:9001` and `-e FILE_LOC=s3://depot/data/grids/hrrr/...`.

### Process local files:

Here's how you'd run a local build of version `v3.0.0` to process a specific forecast hour of HRRR
that you have stored locally:

```shell
docker run --rm --name column-extractor -v ~/data:/data -e FILE_LOC="/data/input/hrrr.t00z.wrfnatf00.grib2 /data/input/hrrr.t00z.wrfprsf00.grib2" -e ENDPOINT_URL=None -e LOG_DIR=/data/output -e OUTPUT_PATH=None column_extractor:v3.0.0
```

> [!NOTE]
> The `-v` mount is now required, so that local files can be exposed within the running container.
> The mounted path must contain `/input` and `/output` subdirectories. `FILE_LOC` will now be
> relative to `/data/input`.

## Running the container in "interactive" mode

This method overrides the Python command defined by `CMD` in the Dockerfile and opens an
interactive terminal session running a bash shell. This allows you to manually run
the Python command and specify command line arguments as you normally would.

```shell
docker run --rm -it --name column-extractor --entrypoint /bin/bash column_extractor:v3.0.0
```
will launch a terminal session within the container from your terminal.

You can then run the application manually:
```shell
(base) column_extractor@dca26dc9cb2a:~$ python3 src/main.py -c ./conf -f s3://noaa-hrrr-bdp-pds/hrrr.20240701/conus/hrrr.t00z.wrfnatf00.grib2 s3://noaa-hrrr-bdp-pds/hrrr.20240701/conus/hrrr.t00z.wrfprsf00.grib2 -l /data/log
```

> [!NOTE]
> The same mount considerations still apply. Though, in this case, the container won't
> automatically exit when the Python process completes. This gives you some time to
> interrogate the output file within the container, if needed.