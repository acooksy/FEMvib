# FEMvib
# solver for the vibrational Schrodinger equation based on the finite element method
# FEMvib requires LIBMESH, which in turn requires compatible versions of PETSC and SLEPC
# The current version is using petsc-3.23.4, slepc-3.23.2, libmesh-1.8.1, eigen-3.4.0.
# 
# The program suite is currently a patchwork of C++ files (the main engine) and support
#  scripts in Python and Perl.  We will convert the perl to python, but haven't gotten there yet.
#
# The tools/ directory contains scripts and executables of varying degrees of importance:
#
#  * Executed by femvib_run.sh (unless output already exists) *
# - tune: the executable of src/tune.c, scans the results/energyP file to obtain Kriging tuning parameters.
# - tune_run.sh: bash script called by femvib_rub.sh to forward arguments and execute tune.
# - paramsearch.pl: perl script to extract optimal Kriging r value from tune results
# - gmat: the executable of src/gmat.c, solves the Eckart conditions and calculates G matrix values
# - gmatN.pl: perl script to export gmat output into separate files for interpolation.
#
#  * Pre- and post-processing scripts *
# - es_parse.py: script to extract data from Gaussian 16 log files, 
#    creates first-pass versions of energyP and coordinates files.
# - extrpars.py: modules for es_parse.py.
# - node_search.py: crude parser of eigenstate*.dat files to estimate number of nodes along each PES coordinate.

