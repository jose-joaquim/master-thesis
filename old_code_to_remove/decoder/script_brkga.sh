#! /bin/bash

make brkga
path_results="../results/brkga/vrbsp/U_"

max_time=10

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 64;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p ${path_results_final}
        echo "./decoder_brkga ${inst} ${v} ${path_results_final} ${max_time}"
        ./bin/decoder_brkga ${inst} ${v} ${path_results_final} ${max_time}

        # sol_file="solution.txt"
        # ./check ${inst} ${v} ${sol_file}
    done;
done;
