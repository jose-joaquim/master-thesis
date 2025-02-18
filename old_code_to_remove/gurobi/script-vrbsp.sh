#! /bin/bash

lower=1
upper=5
for ((v=lower;v<=upper;v++));
do
    for inst in 8;
    do
        var="vrbsp-ybm"

        opath_results="results/${var}/U_"
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p $path_results_final
        ./mainmodule.py ${inst} ${v} ${path_results_final} 1 ${var}
    done;
done;
