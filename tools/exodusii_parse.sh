#!/bin/bash
# 15-Oct-2024  ALCooksy
# This script converts an ExodusII text output file (AFTER converting to text using ncdump)
#  into x,y,psi or x,y,z,psi format (where psi is the eigenstate output, for us the wavefunction)
#  for display by standard plotting programs such as gnuplot.
# The ExodusII output groups each coordinate and output value in a separate comma-delimited list.
#  So we need to find each chunk of data, read the values, and then print back out in the
#  appropriate format.
# For this purpose, the coordinate and psi values are just strings; no numerical processing required.
# 13-Jul-2025
# Revised to incorporate ncdump command and sorting.
if [ "$#" -lt  2 ] ; then
 echo "usage: exodusii_parse.sh  [Exodus txt filename]  [# dimensions]"
 exit
fi
base=${1%.*}
filein=$base.e
fileout=$base.tmp1
ndim=$2  # we need to know how many dimensions we're working with
if [ -e "$fileout" ] ; then
 echo "output file $fileout exists.  Please remove or rename."
 exit
fi

# looks like we'll do this
echo "converting Exodus II file $filein "
ncdump $filein > $fileout
filein=$fileout
fileout=$base.tmp2
if [ -e "$fileout" ] ; then
 echo "output file $fileout exists.  Please remove or rename."
 exit
fi

# convert ncdump output to XYZ file
echo "parsing file $filein "
# find the line numbers in the input file where the value, x, y, and z arrays begin
startline_v=$(($(awk '/vals_nod_.* =/ {print NR} ' "$filein" | xargs)))
startline_x=$(awk '/coordx =/ {print NR} ' "$filein")
startline_y=$(awk '/coordy =/ {print NR} ' "$filein")
if [ "$ndim" -gt 2 ]; then
 startline_z=$(awk '/coordz =/ {print NR} ' "$filein")
fi
# primitive, but I find the end of each block of data by finding all the blank lines,
# then locating the first blank line after each of the start lines above.
str_blanklines=$(awk '/^$/ { ORS=" "; print NR} ' "$filein")
IFS=' ' read -ra blanklines <<< "$str_blanklines"

# get the values line numbers
i="0"
while [ "${blanklines[$i]}" -lt "$startline_v" ]
do
 i=$((i+1))
 endline_v=$(echo "${blanklines[$i]}" | xargs)
 endline_v=$((endline_v))
done
# get the coordx line numbers
i="0"
while [ "${blanklines[$i]}" -lt "$startline_x" ]
do
 i=$((i+1))
 endline_x=$(echo "${blanklines[$i]}" | xargs)
done
# get the coordy line numbers
i="0"
while [ "${blanklines[$i]}" -lt "$startline_y" ]
do
 i=$((i+1))
 endline_y=$(echo "${blanklines[$i]}" | xargs)
done
if [ "$ndim" -gt 2 ]; then
# get the coordz line numbers
i="0"
while [ "${blanklines[$i]}" -lt "$startline_z" ]
do
 i=$((i+1))
 endline_z="${blanklines[$i]}"
done
fi
startline_v=$((startline_v+1))
endline_v=$((endline_v-1))
endline_x=$((endline_x-1))
endline_y=$((endline_y-1))
# echo "values: $startline_v $endline_v "  # debug
# echo "coordx: $startline_x $endline_x "  # debug
# echo "coordy: $startline_y $endline_y "  # debug

# read in the lines
str_values=$(sed -n "$startline_v,$endline_v p" "$filein" | xargs | tr -d '\n;')
IFS=',' read -ra values <<< "$str_values"
str_coordx=$(sed -n "$startline_x,$endline_x p" "$filein" | xargs | tr -d '\n;')
IFS=', ' read -ra coordx <<< "$str_coordx"
str_coordy=$(sed -n "$startline_y,$endline_y p" "$filein" | xargs | tr -d '\n;')
IFS=', ' read -ra coordy <<< "$str_coordy"
if [ "$ndim" -gt 2 ]; then
str_coordz=$(sed -n "$startline_z,$endline_z p" "$filein" | xargs | tr -d '\n;')
IFS=', ' read -ra coordz <<< "$str_coordz"
fi
# This way of assembling the x,y,z,value arrays begins with two bad values in
# the coord arrays (coordx[0] = "coodx" , coordx[1] = "="), and all arrays end with ";"
# so we subtract 1 from #values[@] and 3 from the others to get the correct number of 
# good array values.  And below we start counting the coords from index 2 instead of index 0.
n_values=$((${#values[@]}))
n_coordx=$((${#coordx[@]}-2))
n_coordy=$((${#coordy[@]}-2))
if [ "$ndim" -gt 2 ]; then
n_coordz=$((${#coordz[@]}-2))
fi
echo "$n_values points detected"

if [ $n_coordx -ne $n_coordy ]; then
 echo "number of x values " $n_coordx " not equal to number of y values " $n_coordy
 echo $str_coordx > coordx.out
 echo $str_coordy > coordy.out
 echo "See files coordx.out and coordy.out"
 exit
fi
if [ $n_values -ne $n_coordx ]; then
 echo "number of psi values " $n_values " not equal to number of points " $n_coordx
 echo $str_coordx > coordx.out
 echo $str_values > values.out
 echo "See files coordx.out and values.out"
 exit
fi
# echo "# values = " $n_values   # debug
# echo "# coordx = " $n_coordx   # debug
# echo "# coordy = " $n_coordy   # debug
echo "point 0:  ${coordx[2]}  ${coordy[2]} ${values[0]}  "
echo "point ${n_values}: ${coordx[$n_coordx+1]}  ${coordy[$n_coordy+1]} ${values[$n_values-1]}  "
touch "$fileout"
i="0"
while [ $i -le $n_coordx ]
do
  if [ "$ndim" -gt 2 ]; then
    echo "${coordx[$i+2]}" "${coordy[$i+2]}" "${coordz[$i+2]}" "${values[$i]}" >> "$fileout"
  else
    echo "${coordx[$i+2]}" "${coordy[$i+2]}" "${values[$i]}" >> "$fileout"
  fi
  echo -ne "point  $i \r"
  i=$((i+1))
done
rm $filein

# sort the output file
echo "Sorting file $filein "
filein=$fileout
fileout=$base.dat
if [ -e "$fileout" ] ; then
 echo "output file $fileout exists.  Please remove or rename."
 exit
fi
if [ "$ndim" -gt 2 ]; then
  sort -g -k1,1 -k2,2 $filein > $fileout
else 
  sort -g -k1,1 -k2,2 -k3,3 $filein > $fileout
fi
rm $filein
echo ""
echo "Output saved in file $fileout"
exit
