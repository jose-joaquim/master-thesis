#! /bin/bash

make vrbsp
# make vrcheck
path_results="results/vrbsp/U_"

max_time=15

lower=1
upper=10
for ((v=lower;v<=upper;v++));
do
    for inst in 64;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p ${path_results_final}
        echo "./vrbsp_vns ${inst} ${v} ${path_results_final} ${max_time}"
        ./vrbsp_vns ${inst} ${v} ${path_results_final} ${max_time}

        # sol_file="solution.txt"
        # ./check ${inst} ${v} ${sol_file}
    done;
done;

