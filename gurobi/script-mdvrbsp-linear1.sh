#! /bin/bash

path_results="results/md-vrbsp-linear1/U_"

lower=1
upper=5
for ((v=lower;v<=upper;v++));
do
    for inst in 128 256 512;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./mdvrbsp_linear1.py ${inst} ${v} ${path_results_final}
    done;
done
