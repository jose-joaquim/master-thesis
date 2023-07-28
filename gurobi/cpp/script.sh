#!/bin/sh

for j in 8 16 32 64 128 256 512 1024 2048;
do
    for i in {1..30};
    do
        ./out ${j} ${i} ../../instances/U_${j}/MD-VRBSP_U_${j}_${i}.txt >> obj
    done;

    mkdir -p U_${j};
    mv sol${j}* U_${j};
    mv obj U_${j};
done;
