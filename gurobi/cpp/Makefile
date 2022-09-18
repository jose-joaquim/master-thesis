GUROBI_HOME=/Library/gurobi952/macos_universal2
GUROBI_INC=$(GUROBI_HOME)/include

all:
	g++ -std=c++11 -o out gurobitmp.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi95 -lgurobi_c++

v2:
	g++ -std=c++11 -o out gurobitmp2.cc -I$(GUROBI_INC) -L$(GUROBI_HOME)/lib -lgurobi95 -lgurobi_c++
