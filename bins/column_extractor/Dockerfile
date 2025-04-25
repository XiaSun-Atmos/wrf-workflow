# using the latest micromamba base image
FROM mambaorg/micromamba

# define code version for output files
ARG ARG_CODE_VERSION="v.DEVELOPMENT"

# create user so that *most* actions not performed as root
ARG NEW_MAMBA_USER=column_extractor
ARG NEW_MAMBA_USER_ID=1000
ARG NEW_MAMBA_USER_GID=1000

# as root, create user specified above
USER root
RUN usermod "--login=${NEW_MAMBA_USER}" "--home=/app" \
        --move-home "-u ${NEW_MAMBA_USER_ID}" "${MAMBA_USER}" && \
    groupmod "--new-name=${NEW_MAMBA_USER}" "-g ${NEW_MAMBA_USER_GID}" "${MAMBA_USER}" && \
    # Update the expected value of MAMBA_USER for the _entrypoint.sh consistency check.
    echo "${NEW_MAMBA_USER}" > "/etc/arg_mamba_user" && \
    :
ENV MAMBA_USER=$NEW_MAMBA_USER

# create data + app directories
RUN mkdir -p /data/input /data/output && \
    chown -R $MAMBA_USER:$MAMBA_USER /data && \
    chown -R $MAMBA_USER:$MAMBA_USER /app

USER $MAMBA_USER

# copy in Python requirements file so that Micromamba can build environment
COPY requirements.yml requirements.yml

# build environment as specified by requirements.yml
RUN micromamba install -y -n base -f requirements.yml && \
    # clean micromamba cache
    micromamba clean -y --all --force-pkgs-dirs && \
    # purge pip cache
    /opt/conda/bin/pip cache purge && \
    # delete miscellaneous eccodes files we don't need for reading GRIB2
    rm -r /opt/conda/share/eccodes/definitions/{bufr,taf,metar,grib1,grib3,mars} && \
    # delete Conda metadata we don't need at runtime
    rm -r /opt/conda/conda-meta && \
    # delete all tests we don't need at runtime
    find /opt/conda/lib/python3.9/ -type d -name tests -prune -exec rm -r {} \;

# /app is working directory
WORKDIR /app

# copy source code in
COPY src src

# copy configuration ~stuff~ in (configuration files, profiling sites, grib mapping, etc)
COPY conf conf

# update ECCODES_DEFINITION_PATH to include our custom definitions
ENV ECCODES_DEFINITION_PATH=/app/conf/custom_definitions:/opt/conda/share/eccodes/definitions/

# set code version environment variable
ENV CODE_VERSION $ARG_CODE_VERSION

# run workflow, accepting file location and endpoint URL as environment variables
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
