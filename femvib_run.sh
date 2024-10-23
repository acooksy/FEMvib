#!/bin/bash
# femvib_run.sh runs (a) tune.c to scan the PES for optimal
#  values of the Kriging shape parameters
#  (b) paramsearch.pl to pull those values from results/paramtune
#  (c) FEMvib to create the FEM mesh and integrate the Schrodinger Eq.
# Requires # dimensions as input for tune.c.
# Requires input FEMvib file "input" with default values for all arguments.
if [ "$#" -eq 0 ]; then
  echo "Usage: femvib_run.sh [ndim]    where ndim = # dimensions = 2 or 3"
  exit 1;
elif [ "$1" -ne 2 ] && [ "$1" -ne 3 ]; then 
  echo "Incorrect # dimensions"
  echo "Usage: femvib_run.sh [ndim]    where ndim = # dimensions = 2 or 3"
  exit 1;
fi
ndim=$1

# Run tune if the opp and optr files don't already exist
if [ ! -e "results/opp" ] || [ ! -e "results/optr" ]; then 
  echo "Running tune to generate Kriging interpolation parameters"
  ./tools/tune_run.sh $ndim    # generates results/paramtune
	# extract optimal Kriging r value from results/paramtune
  ./tools/paramsearch.pl ./results/paramtune 
elif [ -e "results/opp" ] && [ -e "results/optr" ]; then
  echo "Files results/opp and results/optr exist.  Skipping tune."
fi

# Run gmat if the matrixg files don't already exist
if [ $ndim -eq 2 ]; then 
  if [ ! -e "results/matrixg.xx" ] || [ ! -e "results/matrixg.xy" ] || [ ! -e "results/matrixg.yy" ]; then 
    echo "Running gmat to calculate G matrix element surfaces"
    ./tools/gmat > ./gfile
    ./tools/gmatN.pl ./gfile $ndim ./
    rm -f ./gfile
  elif [ -e "results/matrixg.xx" ] && [ -e "results/matrixg.xy" ] && [ -e "results/matrixg.yy" ]; then 
    echo "Files results/matrixg.?? exist.  Skipping gmat."
  fi
elif [ $ndim -eq 3 ]; then 
  if [ ! -e "results/matrixg.xx" ] || [ ! -e "results/matrixg.xy" ] || [ ! -e "results/matrixg.yy" ] ||
     [ ! -e "results/matrixg.zx" ] || [ ! -e "results/matrixg.zy" ] || [ ! -e "results/matrixg.zz" ]; then 
    echo "Running gmat to calculate G matrix element surfaces"
    ./tools/gmat > ./gfile
    ./tools/gmatN.pl ./gfile $ndim ./
    rm -f ./gfile
  elif [ -e "results/matrixg.xx" ] && [ -e "results/matrixg.xy" ] && [ -e "results/matrixg.yy" ] &&
       [ -e "results/matrixg.zx" ] && [ -e "results/matrixg.zy" ] && [ -e "results/matrixg.zz" ]; then 
    echo "Files results/matrixg.?? exist.  Skipping gmat."
  fi
fi

echo "Running FEMvib"
# If you get odd termination of FEMvib, first check that input has correct format.
#  Must identify values for *all* arguments; elemTypeStr must match one of the
#  available FEM types for that dimensionality.
./femvib -f results/input 
# clean up temporary files
# rm -f results/ndim results/opp results/optr
exit 0 ;

