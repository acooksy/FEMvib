#!/usr/bin/env bash
 dim=$1
# np='10 20 30'
# rpar='0.6 1.0 1.4 1.8 2.2 2.6 3.0'
 np='20 25 30 35 40 45 50 60 70 80 90 100'
 rpar='1.1 1.5 2.1 2.5 3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 8.0'

if [ -e ./results/paramtune ]; then
       	rm -f ./results/paramtune
fi
for j in $np
 do
 for i in $rpar
        do
	echo -ne "i,j = $i $j \033[0K\r"
        ./tools/tune $dim $i $j >> ./results/paramtune
        done
 done
echo ""

