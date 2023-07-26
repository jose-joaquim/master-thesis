vns:
	g++ -Wall -std=c++11 VNS.h

vrbsp:
	g++ -Wall -std=c++2a -g -o vrbsp_vns vrbsp_vns.cpp

mdvrbsp:
	mkdir -p bin
	g++ -Wall -DMDVRBSP -O0 -g -std=c++2a mdvrbsp.cpp -o bin/out
