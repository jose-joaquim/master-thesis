#! /bin/bash

make brkga

for inst in 64 128 512 1024 256 32 16 8; do
	for v in 25 7 28 19 1 5 29 15 21 30 14 9 16 23 24 11; do
		./bin/decoder_brkga ../../instances/dr/U_${inst}/MD-VRBSP_U_${inst}_${v}.txt
	done
	echo "done ${inst}"
done
