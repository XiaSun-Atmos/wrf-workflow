###
### This script is included in /src/tests. We can call it when running
### the container as an alternative to running the actual processing code.
### By doing so, it runs the unit tests included in /src/tests.
### This is done in practice via /integration/test_build.sh.
###

# change into src/ directory
cd src/

# run all tests under src/tests
/opt/conda/bin/python3 -m unittest -v
