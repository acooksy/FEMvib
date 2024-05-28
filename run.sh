#!/bin/bash

#set -x

if [ -f "make.log" ] ; then
	rm -f make.log
fi

## source $LIBMESH_DIR/examples/run_common.sh
source ./run_common_.sh

## example_name=eigenproblems_ex2
example_name=femvib_ep-ex2

## options="-n 5 -eps_type lapack"
## options="-n 10 -nx 10 -ny 10 -nz 10"
## options="-ndim 3 -n 10 -nx 10 -ny 10 -nz 10 -x0 0.5 -y0 0.5 -z0 0.5  -elemTypeStr HEX27"  ## 3D
options="-ndim 2 -n 10 -nx 40 -ny 40 -x0 1 -y0 1 -elemTypeStr QUAD9"  ## 2D
## options="-ndim 1 -n 10 -nx 40 -x0 1 -elemTypeStr EDGE2"  ## 1D

run_example "$example_name" "$options"

# No benchmark_example here - it spends 99.8+% of its time in the
# SLEPc solve at even small scale, with that eps_type
