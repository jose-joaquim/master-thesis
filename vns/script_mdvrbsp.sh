#! /bin/bash

make mdvrbsp
# make mdcheck
path_results="results/vns-vrbsp/U_"

max_time=10

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 16;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p ${path_results_final}
        echo "./mdvrbsp_vns ${inst} ${v} ${path_results_final} ${max_time}"
        ./mdvrbsp_vns ${inst} ${v} ${path_results_final} ${max_time}

        # sol_file="solution.txt"
        # ./check ${inst} ${v} ${sol_file}
    done;
done;
