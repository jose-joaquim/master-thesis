GUROBI_HOME=/Library/gurobi1002/macos_universal2
GUROBI_INC=$(GUROBI_HOME)/include

all:
	mkdir -p bin
	g++ -std=c++2a -DMDVRBSP -o bin/out mdvrbsp.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi100 -lgurobi_c++

vrbsp:
	mkdir -p bin
	g++ -std=c++2a -g -o bin/out vrbsp.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi100 -lgurobi_c++
