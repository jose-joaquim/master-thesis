path_results="results/md-vrbsp/U_"
for inst in 8;
do
    lower=1
    upper=2
    for ((v=lower;v<=upper;v++));
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./mdvrbsp_linear1.py ${inst} ${v} ${path_results_final}
    done;
done
