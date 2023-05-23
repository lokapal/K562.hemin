#!/bin/bash
#filename=$(basename "$1")
#extension="${filename##*.}"
#filen="${filename%.*}"
#outname=$filen
#outname+=".bedGraph"
bedops --partition $1 | bedmap --echo --echo-map-id-uniq --delim '\t' - $1 | awk '{ n = split($4, a, ";"); sum = 0; for(i = 1; i <= n; i++) { sum += a[i]; } meann = sum/n; print $1"\t"$2"\t"$3"\t"meann; }' - > $2
