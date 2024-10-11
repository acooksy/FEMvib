#!/bin/bash
ndim=$1
neig=$2
nx=0
ny=0
nz=0
x0=0
x1=0
y0=0
y1=0
z0=0
z1=0
elemtyp=""
if [[ $ndim == 1 ]] ; then
	nx=$3
	x0=$4
	x1=$5
	elemtyp=$6 ;
elif [[ $ndim == 2 ]] ; then
	nx=$3
	ny=$4
	x0=$5
	x1=$6
	y0=$7
	y1=$8
	elemtyp=$9 ;
elif [[ $ndim == 3 ]] ; then
	nx=$3
	ny=$4
	nz=$5
	x0=$6
	x1=$7
	y0=$8
	y1=$9
	z0=${10}
	z1=${11}
	elemtyp=${12} ;
fi

echo "Running gmat"
./tools/gmat > ./gfile
./tools/gmatN.pl ./gfile $ndim ./
rm -f ./gfile
echo "FEMvib"
./femvib -ndim $ndim -n $neig -nx $nx -ny $ny -nz $nz -x0 $x0 -x1 $x1 -y0 $y0 -y1 $y1 -z0 $z0 -z1 $z1 -elemTypeStr $elemtyp
exit 0 ;

