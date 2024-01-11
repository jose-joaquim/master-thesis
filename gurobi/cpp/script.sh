#!/bin/sh

for j in 64; 
do
    for i in {1..1};
    do
        ./bin/out ${j} ${i} ../../instances/U_${j}/MD-VRBSP_U_${j}_${i}.txt gurobi_sol${j}_${i}.sol
    done;

    # mkdir -p U_${j};
    # mv sol${j}* U_${j};
    # mv obj U_${j};
done;
