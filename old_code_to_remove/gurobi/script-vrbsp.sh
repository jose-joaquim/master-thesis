#! /bin/bash

lower=1
upper=1
for ((v = lower; v <= upper; v++)); do
	for inst in 8; do
		var="vrbsp-ybm"

		opath_results="results/${var}/U_"
		path_results_final=${path_results}${inst}
		echo $path_results_final
		mkdir -p $path_results_final
		./mainmodule.py ${path} ${var}
	done
done
