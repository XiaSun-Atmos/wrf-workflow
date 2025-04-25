###
### This script is included in /integration. It takes a "truth" build (via tag)
### and a "test" build (via tag) as input. It'll then fetch GFS/RAP/HRRR test files
### and process them using the truth and test builds. It'll then run our unit
### tests against the output files to verify consistency between builds.
###
### To compare v1.1.0 (development) against v1.0.3 (truth), run:
###
### "bash ./integration/test_build.sh -t v1.0.3 -d v1.1.0"
###
### (include optional "-l yes" if you have local files to test)
###

echo "Testing column extraction workflow..."
echo ""

### get command line arguments
while getopts t:d:l: flag

do

    case "${flag}" in
        t) truth_build=${OPTARG};; # the container image (tag) denoting the "truth" build
        d) dev_build=${OPTARG};;   # the container image (tag) denoting the "development" build
        l) local=${OPTARG};;       # whether or not to fetch test files from remote source
    esac

done

echo "Development build: $dev_build"
echo "Truth build:       $truth_build"
echo ""

### step 1: create necessary directories if don't exist

# make input/output/testing directories
mkdir -p ./development/data/{input,output,testing}

# declare models, buckets, paths, objects to pull
models=("rap" "hrrr" "gfs")
buckets=("noaa-rap-pds" "noaa-hrrr-bdp-pds" "noaa-gfs-bdp-pds")
paths=("rap.20240701/" "hrrr.20240701/conus/" "gfs.20240701/00/atmos/")
objects=(
    "--include rap.t00z.wrfnatf00.grib2 --include rap.t00z.wrfnatf03.grib2"
    "--include hrrr.t00z.wrfnatf00.grib2 --include hrrr.t00z.wrfprsf00.grib2 --include hrrr.t00z.wrfnatf03.grib2 --include hrrr.t00z.wrfprsf03.grib2"
    "--include gfs.t00z.pgrb2.0p25.f000 --include gfs.t00z.pgrb2b.0p25.f000 --include gfs.t00z.pgrb2.0p25.f003 --include gfs.t00z.pgrb2b.0p25.f003"
)

for idx in "${!models[@]}"

do

    echo "Testing ${models[$idx]}..."
    echo ""

    ### step 2: fetch model files from /public or use local files
    ### (this includes an FH00 and FHNN file, since behavior can be different)

    if ! [ $local ]; then

        echo "Fetching ${models[$idx]} GRIB output to test against..."
        echo ""

        aws s3 cp s3://${buckets[$idx]}/${paths[$idx]} ./development/data/input/${models[$idx]}/ --recursive --exclude "*" ${objects[$idx]} --no-sign-request --no-progress
        echo ""

    fi

    for forecast_hour in "00" "03"

    do

        files=$(find ./development/data/input/${models[$idx]} -regex ".*f0*${forecast_hour}$" -o -regex ".*f0*${forecast_hour}.grib2" | sed "s#./development##g")

        ### step 3: run build against which you want to test (aka "truth")

        echo "Running $truth_build build on $files..."
        echo ""

        # run container
        docker run --rm --name column-extractor -v ./development/data:/data -e FILE_LOC="${files}" -e ENDPOINT_URL=None -e LOG_DIR=/data/output -e OUTPUT_PATH=None column_extractor:$truth_build

        # move output file to ./development/data/testing as truth_output.nc
        mv ./development/data/output/profiler_sites* ./development/data/testing/truth_output.nc

        ### step 4: run development build

        echo "Running $dev_build build on $files..."
        echo ""

        # run container
        docker run --rm --name column-extractor -v ./development/data:/data -e FILE_LOC="${files}" -e ENDPOINT_URL=None -e LOG_DIR=/data/output -e OUTPUT_PATH=None column_extractor:$dev_build

        # move output file to ./development/data/testing as truth_output.nc
        mv ./development/data/output/profiler_sites* ./development/data/testing/test_output.nc

        ### step 5: run tests

        echo "Running tests..."
        echo ""

        docker run --rm --name column-extractor -v ./development/data:/data column_extractor:$truth_build ./src/tests/run_tests.sh
        echo ""

    done

done
