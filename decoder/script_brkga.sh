#! /bin/bash

make
path_results="results/vrbsp/U_"

max_time=20

lower=1
upper=2
for ((v=lower;v<=upper;v++));
do
    for inst in 64;
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        mkdir -p ${path_results_final}
        echo "./decoder_brkga ${inst} ${v} ${path_results_final} ${max_time}"
        ./decoder_brkga ${inst} ${v} ${path_results_final} ${max_time}

        # sol_file="solution.txt"
        # ./check ${inst} ${v} ${sol_file}
    done;
done;
