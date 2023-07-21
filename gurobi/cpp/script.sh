#!/bin/sh

U=64

for i in {29..29};
do
  ./out ${U} ${i} ../../../instances/U_${U}/MD-VRBSP_U_${U}_${i}.txt
done
