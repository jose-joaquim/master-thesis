#! /bin/bash

path_results="results/vrbsp-mauricio/U_"

lower=1
upper=30
for ((v=lower;v<=upper;v++));
do
    for inst in 512 1024 2048;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./vrbsp_mauricio.py ${inst} ${v} ${path_results_final}
    done;
done
