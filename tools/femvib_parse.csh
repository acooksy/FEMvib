#!/bin/csh -f
# 18-Jun-2014
# This script reads a gaussian DFT, HF, or CC log file and extracts optimized geometries 
#    and corresponding energies.
# If it detects "Stationary point" anywhere in the file, then it assumes that all the desired geometries 
#    are fully optimized (i.e., it should ignore any incomplete optimizations at the end of the file).
# If it detects no "Stationary point", then it assumes that only the last energy and geometry desired.
# If an optional second filename is included, this is assumed to be a single geometry optimization job;
#  all energies are then calculated relative to the energy of the optimized geometry and the optimized
#  geometry is placed at the head of the coordinates file as the reference geometry.
# If there is no optimization information, the absolute energies in Eh are printed instead.
#
#  revisions:
#  31-Jul-2014 If no stationary points are reported, script will instead extract the last energy and 
#    corresponding geometry.  
#    Also added auto-detect number of atoms.
#    Fixed bug in values array.
#    Energies files also includes rotational constants.
#  15-Sep-2024
#    Add optional second argument for filename of opt job to set energy minimum and reference geometry
#    Convert energies to cm-1 in energies file
#    Add ref geometry to geoms file

# Check that we have one argument; else print syntax
if (( $#argv != 1) && ( $#argv != 2 )) then
 echo "syntax: femvib_parse.csh [scan output filename] {opt output filename}"
 exit(1)
endif

# Initialize variables
set file = $1
set flag_opt = 0
if ( $#argv > 1) then
 set flag_opt = 1
 set optfile = $2
 if ( ! `grep -c "Stationary" $optfile `) then
  echo "No optimized geometry found in ${optfile}.  Exiting."
  exit(1)
 endif
endif
set flag_cc = 0
if (`grep -c "E(CORR)" $file`) set flag_cc = 1
set flag_zm = 1
if (`grep -c "Name  Definition" $file`) set flag_zm = 0
set n_atoms = `grep NAtoms $file | tail -1 | awk '{print $2}'`

# The next 3 lines assemble arrays of the line numbers in $file that contain the strings
#  "Station" (indicates an optimized geometry found), "E(" (for the energy), "Standard" (for the geometry)
#  For a CC or QCI run, instead search for "E(CORR)" for the energy (otherwise may get E(Z))
set lines_stat = `awk '/Station/ {print NR} ' $file`
if ( $flag_opt ) set optlines_stat = `awk '/Station/ {print NR} ' $optfile`
if ( $flag_cc ) then
 echo "coupled cluster job"
 set lines_energy = `awk '/E\(CORR)/ {print NR} ' $file`
 if ( $flag_opt ) set optlines_energy = `awk '/E\(CORR)/ {print NR} ' $optfile`
else if ( ! $flag_cc ) then
 set lines_energy = `awk '/E\(/ {print NR} ' $file`
 if ( $flag_opt ) set optlines_energy = `awk '/E\(/ {print NR} ' $optfile`
endif
set lines_stand = `awk '/Standard orientation/ {print NR} ' $file`
if ( $flag_opt ) set optlines_stand = `awk '/Standard orientation/ {print NR} ' $optfile`
set lines_rot = `awk '/Rotational constants/ {print NR} ' $file`
if ( $flag_opt ) set optlines_rot = `awk '/Rotational constants/ {print NR} ' $optfile`
set num_stat = $#lines_stat
set num_energy = $#lines_energy
set num_stand = $#lines_stand
set num_rot = $#lines_rot
set i=1
set lines_E_opt = 0
set lines_r_opt = 0
while ($i <= $num_stat)
 set lines_E_opt = ($lines_E_opt 0)
 set lines_r_opt = ($lines_r_opt 0)
 @ i++
end
set i=1
set lines_g_opt = 0
while ($i < $num_stat)
 set lines_g_opt = ($lines_g_opt 0)
 @ i++
end

if ( $num_energy == 0 ) then
 echo "No energy values found in file ${file}.  Exiting."
 exit(1)
endif

echo "number of atoms = $n_atoms" 
echo "number of stationary points = $num_stat" 
echo "number of energies = $num_energy" 
echo "number of geometries = $num_stand" 

# Check that the output file doesn't already exist
set flag_exit = 0
if (-e ${file}.energies ) then
 echo "Output file ${file}.energies already exists.  Delete or rename, then restart."
 set flag_exit = 1
endif
if (-e ${file}.geoms) then
 echo "Output file ${file}.geoms already exists.  Delete or rename, then restart."
 set flag_exit = 1
endif
if ( $flag_exit ) then 
 echo "Exiting."
 exit(2)
endif

# Identify the PES coordinates by finding parameters tagged "Frozen" or "Scan".
#  params[i] i=1..num_params
set params = `awk '/   Scan   / || /Frozen/ || /frozen, / {print $2}' $file`
set num_params = $#params
echo "Number of PES coordinates detected = $num_params"
if ($num_params == 0) then
 echo "No coordinates found.  If this is a restart job, that may be the problem.  Exiting."
 exit(1)
endif
echo "PES coordinates: $params "
# Check if any of these are angular coordinates ; in that case we have to convert degrees to radians
set i=1
set flag_ang = ()
while ($i <= $num_params)
 set flag_ang = ($flag_ang 0)
 set char = `echo $params[$i] | grep -o '^.'`
 if ( $char == "D") set flag_ang[$i] = 1
 if ( $char == "A") set flag_ang[$i] = 1
 @ i++
end

# Find the values of the frozen PES coordinates
# First initialize arrays
set i=1
set values = ()
while ($i <= $num_params)
 set values = ($values 0)
 if ($flag_opt) set optvalues = ($values 0)
 @ i++
end
# For stationary points, pull values from the Optimized Geometry output
if ($num_stat > 0 ) then 
 set i = 1
 while ($i <= $num_params )
  if ( ! $flag_zm ) then 
   if ( ! $flag_ang[$i] ) then
    set values[$i] = `grep "! .*$params[$i] .* -DE/DX" $file | awk ' {print $4}' `
    if ($flag_opt) set optvalues[$i] = `grep "! .*$params[$i] .* -DE/DX" $optfile | awk ' {print $4}' `
   else if ( $flag_ang[$i] ) then 
    set values[$i] = `grep "! .*$params[$i] .* -DE/DX" $file | awk '{ val = $4 * 0.0174533 ; print val }' `
    if ($flag_opt) set optvalues[$i] = `grep "! .*$params[$i] .* -DE/DX" $optfile | awk '{ val = $4 * 0.0174533 ; print val}' `
   endif
  else if ( $flag_zm ) then 
   if ( ! $flag_ang[$i] ) then
    set values[$i] = `grep "! .*$params[$i] .* -DE/DX" $file | awk ' {print $3}' `
    if ($flag_opt) set optvalues[$i] = `grep "! .*$params[$i] .* -DE/DX" $optfile | awk ' {print $3}' `
   else if ( $flag_ang[$i] ) then 
    set values[$i] = `grep "! .*$params[$i] .* -DE/DX" $file | awk '{ val = $3 * 0.0174533 ; print val}' `
    if ($flag_opt) set optvalues[$i] = `grep "! .*$params[$i] .* -DE/DX" $optfile | awk '{ val = $3 * 0.0174533 ; print val}' `
   endif
  endif
  @ i++
 end
# If there's no stationary point, then we have to pull the values from the top of the log file
else if ($num_stat == 0 ) then
 set i = 1
 while ($i <= $num_params )
# need to differentiate opt=z-matrix jobs which have different format here
  if ( ! $flag_zm ) then 
   set values[$i] = `grep "! .*$params[$i] .*rozen\|! .*$params[$i] .* Scan" $file | awk ' {print $4}' `
  else if ( $flag_zm ) then 
   set values[$i] = `grep "! .*$params[$i] .*rozen" $file | awk ' {print $3}' `
  endif
  @ i++
 end
endif
set i = 1
while ($i <= $num_params )
 echo "values[$i] = $values[$i] "
 @ i++
end

# This loop finds all the lines reporting the energy of the optimized geometry
if ($num_stat > 0) then
 set i = 1
 set j = 1
 while (($j <= $num_stat) && ($i <= $num_energy))
  if ( $lines_energy[$i] >= $lines_stat[$j] ) then
   @ j++
  else 
   set lines_E_opt[$j] = $lines_energy[$i]
   @ i++
  endif
 end
endif

# Find lines with energies and rotational constants in optimziation file if any
if ($flag_opt) then
 set optline_E_opt = $optlines_energy[$#optlines_energy]
 set optline_r_opt = $optlines_rot[$#optlines_rot]
endif

# This loop finds all the initial lines of the optimized geometries.
echo "Finding initial lines of optimized geometries"
if ($num_stat > 0) then
set i = 1
set j = 1
 while (($j <= $num_stat) && ($i <= $num_stand))
  if ( $lines_stand[$i] >= $lines_stat[$j] ) then
   @ j++
  else 
   set lines_g_opt[$j] = $lines_stand[$i]
   @ lines_g_opt[$j] += 5
   @ i++
  endif
 end
else if ($num_stat == 0) then
#  we want the index in the next line to be $num_energy rather than $num_stand, 
#   because there may be an extra geometry printed after the energy
 set lines_g_opt[1] = $lines_stand[$num_energy]
 @ lines_g_opt[1] += 5
endif

# Repeat for the optimization file if any
if ($flag_opt) then
 set optline_g_opt = $optlines_stand[$#optlines_stand]
 @ optline_g_opt += 5
endif
echo "Done finding initial lines of optimized geometries"

# Now loop over the results and print the energies and geometries.
echo "Printing energies and geometries"
if ($num_stat == 0) set num_geom = 1
if ($num_stat > 0) set num_geom = $num_stat
touch ${file}.geoms
touch ${file}.energies
set num_geom_tot = $num_geom
if ($flag_opt) @ num_geom_tot ++
echo $num_params $n_atoms $num_geom_tot >> ${file}.geoms
echo "coordinates =  $params" >> ${file}.geoms
echo "coordinates =  $params" >> ${file}.energies

# First parse the optimization file if any
if ($flag_opt) then
 set values_list = () 
 set i = 1
 while ($i <= $num_params)
  set values_list = ($values_list `echo $optvalues[$i] `) 
  @ i++
 end
 set optline_g_end = $optline_g_opt
 @ optline_g_end += $n_atoms
# first print the optimized geom as the reference geometry, includes atom specification
 awk -v num=$optline_g_opt -v num2=$optline_g_end 'NR >= num && NR < num2 {print $2," ", $4," ", $5," ", $6}' $optfile >> ${file}.geoms
# now print optimized geom again as the first of the PES geom list; no atom spec now
 echo $values_list >> ${file}.geoms
 awk -v num=$optline_g_opt -v num2=$optline_g_end 'NR >= num && NR < num2 {print $4," ", $5," ", $6}' $optfile >> ${file}.geoms
 if ( ! $flag_cc ) set optenergy = `awk -v num=$optline_E_opt 'NR == num {print $5}' $optfile`
 if (   $flag_cc ) set optenergy = `awk -v num=$optline_E_opt 'NR == num {print $4}' $optfile`
 set rot = `awk -v num=$optline_r_opt 'NR == num {print $4 " " $5 " " $6}' $optfile`
 set energy = 0.0
# printing rotational constants is an unnecessary option
# set values_list = ( $values_list $energy $rot)
 set values_list = ( $values_list $energy )
 echo $values_list >> ${file}.energies
echo "Printed optimized geometry and energy"
endif

# Now parse the scan file
echo "Parsing scan file"
set j = 1
while ($j <= $num_geom)
 set values_list = () 
 set i = 1
 while ($i <= $num_params)
  set values_list = ($values_list `echo $values[$i] | awk -v num=$j '{print $num}'`) 
  @ i++
 end
 echo $values_list >> ${file}.geoms
 if ( $flag_opt) then 
#  if ( ! $flag_cc ) set energy = `awk -v num=$lines_E_opt[$j] -v optenergy=$optenergy 'NR == num { val = ($5 - optenergy)*219475 ; print val}' $file`
  if ( ! $flag_cc ) then
    set energy = `awk -v num=$lines_E_opt[$j] -v optenergy=$optenergy 'NR == num { val = ($5 - optenergy)*219475 ; print val}' $file`
  endif
  if (   $flag_cc ) set energy = `awk -v num=$lines_E_opt[$j] -v optenergy=$optenergy 'NR == num { val = ($4 - optenergy)*219475 ; print val}' $file`

 else if ( ! $flag_opt) then
  if ( ! $flag_cc ) set energy = `awk -v num=$lines_E_opt[$j] 'NR == num { print $5 }' $file`
  if (   $flag_cc ) set energy = `awk -v num=$lines_E_opt[$j] 'NR == num { print $4 }' $file`
 endif 
 set rot = `awk -v num=$lines_r_opt[$j] 'NR == num {print $4 " " $5 " " $6}' $file`
# printing rotational constants is an unnecessary option
# set values_list = ( $values_list $energy $rot)
 set values_list = ( $values_list $energy )
 echo $values_list >> ${file}.energies
 set lines_g_end = $lines_g_opt[$j]
 @ lines_g_end += $n_atoms
# Correct format for FEMvib coordinates files does not include atom spec except for the reference geometry
# awk -v num=$lines_g_opt[$j] -v num2=$lines_g_end 'NR >= num && NR < num2 {print $2," ", $4," ", $5," ", $6}' $file >> ${file}.geoms
 awk -v num=$lines_g_opt[$j] -v num2=$lines_g_end 'NR >= num && NR < num2 {print $4," ", $5," ", $6}' $file >> ${file}.geoms
 @ j++
 echo -n "geometry " $j " of " $num_geom
end

exit(0)
