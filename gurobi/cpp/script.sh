#!/bin/sh

U=64

for i in {1..1};
do
  ./bin/out ${U} ${i} ../../instances/U_${U}/MD-VRBSP_U_${U}_${i}.txt
done
