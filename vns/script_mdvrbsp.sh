#! /bin/bash

make mdvrbsp
# make mdcheck
path_results="results/mdvrbsp/U_"

max_time=600

lower=9
upper=9
for ((v=lower;v<=upper;v++));
do
    for inst in 512;
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
