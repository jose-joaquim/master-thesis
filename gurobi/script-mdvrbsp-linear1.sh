#! /bin/bash

path_results="results/md-vrbsp-linear1/U_"

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 8;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./mdvrbsp_linear1.py ${inst} ${v} ${path_results_final} 0
    done;
done
