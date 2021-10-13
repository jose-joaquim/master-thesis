#! /bin/bash



lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 8;
    do
        path_results="results/vrbsp/U_"
	    path_results_final=${path_results}${inst};
	    echo $path_results_final
	    ./vrbsp_lin_impr.py ${inst} ${v} ${path_results_final}
    done;

    for inst in 512 1024 2048;
    do
        path_results="results/vrbsp-mauricio/U_"
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./vrbsp_mauricio.py ${inst} ${v} ${path_results_final}
    done;
done;
