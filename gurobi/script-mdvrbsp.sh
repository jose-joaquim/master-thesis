#! /bin/bash

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    # for inst in 16;
    # do
    #     path_results="results/mdvrbsp-quad/U_"
    #     path_results_final=${path_results}${inst};
    #     echo $path_results_final
    #     mkdir -p $path_results_final
    #     ./mainmodule.py ${inst} ${v} ${path_results_final} 1 mdvrbsp-quad
    # done;

    for inst in 64;
    do
        var="mdvrbsp-quad"
        path_results="results/${var}/U_"
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p $path_results_final
        ./mainmodule.py ${inst} ${v} ${path_results_final} 1 ${var}
    done;
done
