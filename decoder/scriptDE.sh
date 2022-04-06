#! /bin/bash

make de
path_results="results/vrbsp/U_"

max_time=20

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 128;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p ${path_results_final}
        echo "./decoder_de ${inst} ${v} ${path_results_final} ${max_time}"
        ./decoder_de ${inst} ${v} ${path_results_final} ${max_time}
    done;
done;
