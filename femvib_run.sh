#!/bin/bash
# femvib_run.sh runs (a) tune.c to scan the PES for optimal
#  values of the Kriging shape parameters
#  (b) paramsearch.pl to pull those values from results/paramtune
#  (c) FEMvib to create the FEM mesh and integrate the Schrodinger Eq.
# Input files should all be in ./results/ ; this is also where output will be printed.
# Requires # dimensions as input for tune.c.
# Requires input FEMvib file "input" with default values for all arguments.

# Read the number of dimensions from the results/input file
ndim=`grep ndim results/input | awk -F '=' '{print $2}'`
if [ "$ndim" -ne 1 ] && [ "$ndim" -ne 2 ] && [ "$ndim" -ne 3 ]; then 
  echo "Incorrect # dimensions"
  echo "Usage: femvib_run.sh [ndim]    where ndim = # dimensions = 1, 2, or 3"
  exit 1;
fi
echo "This is a ${ndim}-dimensional case."

# Run tune if the opp and optr files don't already exist
if [ "$ndim" -ne 1 ]; then
  if [ ! -e "results/opp" ] || [ ! -e "results/optr" ]; then 
    echo "Running tune to generate Kriging interpolation parameters"
    ./tools/tune_run.sh $ndim    # generates results/paramtune
  	# extract optimal Kriging r value from results/paramtune
    ./tools/paramsearch.pl ./results/paramtune 
    echo ""
  elif [ -e "results/opp" ] && [ -e "results/optr" ]; then
    echo "Files results/opp and results/optr exist.  Skipping tune."
  fi
fi

# Run gmat if the matrixg files don't already exist
if [ "$ndim" -eq 1 ]; then 
  if [ ! -e "results/matrixg.xx" ]; then 
    echo "Running gmat to calculate G matrix element surface"
    ./tools/gmat > ./gfile
    ./tools/gmatN.pl ./gfile $ndim ./
    rm -f ./gfile
  elif [ -e "results/matrixg.xx" ]; then 
    echo "File results/matrixg.xx exists.  Skipping gmat."
  fi
elif [ "$ndim" -eq 2 ]; then 
  if [ ! -e "results/matrixg.xx" ] || [ ! -e "results/matrixg.xy" ] || [ ! -e "results/matrixg.yy" ]; then 
    echo "Running gmat to calculate G matrix element surfaces"
    ./tools/gmat > ./gfile
    ./tools/gmatN.pl ./gfile $ndim ./
    rm -f ./gfile
  elif [ -e "results/matrixg.xx" ] && [ -e "results/matrixg.xy" ] && [ -e "results/matrixg.yy" ]; then 
    echo "Files results/matrixg.* exist.  Skipping gmat."
  fi
elif [ "$ndim" -eq 3 ]; then 
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
# If you get odd termination of FEMvib, first check that the file "input" has correct format.
#  One must identify values for *all* arguments, and elemTypeStr must match one of the
#  available FEM types for that dimensionality.
./femvib -f ./results/input 
# clean up temporary files; in general we want to keep opp and optr because "tune" can
# take a long time to run in 3D and/or if there are a lot of grid points.
rm -f ./results/ndim  ./results/paramtune 

exit 0 ;

