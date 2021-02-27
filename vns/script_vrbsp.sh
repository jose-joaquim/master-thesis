#! /bin/bash

make vrbsp
# make vrcheck
path_results="results/vns-vrbsp/U_"

max_time=10

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 8;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p ${path_results_final}
        ./vrbsp_vns ${inst} ${v} ${path_results_final} ${max_time}

        # sol_file="solution.txt"
        # ./check ${inst} ${v} ${sol_file}
    done;
done;

