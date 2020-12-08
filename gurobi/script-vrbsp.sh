path_results="results/vrbsp/U_"
for inst in 8;
do
    lower=1
    upper=1
    for ((v=lower;v<=upper;v++));
    do
        path_results_final=${path_results}${inst};
        echo $path_results_final
        ./vrbsp_new.py ${inst} ${v} ${path_results_final}
    done;
done
