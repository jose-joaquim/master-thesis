#! /bin/bash

make mdvrbsp

arr_of=("MINMAX") # "SUMVIO")
arr_ch=("RANDOM") # "GREEDY")

path_results="results/mdvrbsp/U_"

max_time=60

lower=1
upper=1
for ((v=lower;v<=upper;v++));
do
    for inst in 128;
    do
        for of in "${arr_of[@]}";
        do
            for ch in "${arr_ch[@]}";
            do
                # echo "$of $ch"
                path_var="results/mdvrbsp_${of}_${ch}/U_$inst"
                
                echo $path_var
                mkdir -p $path_var

                echo "./mdvrbsp_vns ${inst} ${v} ${path_var} ${max_time} $of $ch"
                ./mdvrbsp_vns ${inst} ${v} ${path_var} ${max_time} ${of} ${ch}
                # ./mdvrbsp_vns ${inst} ${v}
            done;
        done;

        
        # path_results_final=${path_results}${inst};
        # echo $path_results_final
        # mkdir -p ${path_results_final}
        # 
        # # echo "./mdvrbsp_vns ${inst} ${v} ${path_results_final} ${max_time}"
        # # ./mdvrbsp_vns ${inst} ${v} ${path_results_final} ${max_time} params.prm
        # 
        # ./mdvrbsp_vns ${inst} ${v}
    done;
done;
