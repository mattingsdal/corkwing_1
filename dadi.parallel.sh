#!/bin/bash

for i in `seq 1 1 25`; do

# test models demo dadi #

python2.6 dadi.run.py -f ../Data/SFS.fs -y South -x West -p 50,60,70 -m SI,IM,AM,SC,DB -l -v

((i++))
done

exit 0;

