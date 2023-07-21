GUROBI_HOME=/Library/gurobi1002/macos_universal2
GUROBI_INC=$(GUROBI_HOME)/include

all:
	g++ -std=c++2a -o out gurobitmp.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi100 -lgurobi_c++

v2:
	g++ -std=c++2a -o out gurobitmp2.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi100 -lgurobi_c++

vrbsp:
	mkdir bin
	g++ -std=c++2a -g -o bin/out vrbsp.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi100 -lgurobi_c++
