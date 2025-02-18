GUROBI_HOME=/Library/gurobi1101/macos_universal2
GUROBI_INC=$(GUROBI_HOME)/include

all:
	mkdir -p bin
	g++ -g -std=c++2a -o bin/out main.cc common.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi110 -lgurobi_c++
