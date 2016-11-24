#!/bin/bash

for i in `seq 1 1 25`; do

# test models demo dadi #

python2.6 script_inference_anneal3_newton.py -f ../Data/scandinavia_LD_project40_spectrum_from_dadi.v6.fs -y South -x West -p 50,60,70 -m SI,IM,AM,SC,DB -l -v

((i++))
done

exit 0;

