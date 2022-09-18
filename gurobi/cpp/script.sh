#!/bin/sh

U=8

for i in {1..30};
do
  ./out ${U} ${i} < ../../instances/U_${U}/MD-VRBSP_U_${U}_${i}.txt > /dev/null
done
