#! /bin/bash

path_results="results/md-vrbsp-linear2/U_"

lower=1
upper=30
for ((v=lower;v<=upper;v++));
do
    for inst in 8 16 32 64;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./mdvrbsp_linear2.py ${inst} ${v} ${path_results_final}
    done;
done
