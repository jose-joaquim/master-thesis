GUROBI_HOME=/Library/gurobi1101/macos_universal2
GUROBI_INC=$(GUROBI_HOME)/include

all:
	mkdir -p bin
	g++ -g -std=c++2a -DMDVRBSP -o bin/out mdvrbsp_columngeneration.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi110 -lgurobi_c++

vrbsp:
	mkdir -p bin
	g++ -std=c++2a -DORIG_OF -g -o bin/out vrbsp.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi110 -lgurobi_c++
