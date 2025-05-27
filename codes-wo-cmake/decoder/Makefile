OMP_PREFIX="/opt/homebrew/Cellar/libomp/20.1.3"
JSON_PATH="./include/single_include"

brkga: decoder.o
	g++ -std=c++20 -c -g -L${OMP_PREFIX}/lib -I${OMP_PREFIX}/include -I../../include/ -I${JSON_PATH} HeuristicBRKGA.cpp
	g++ -std=c++20 -o ./bin/decoder_brkga HeuristicBRKGA.o basic.o

decoder.o:
	mkdir -p bin
	g++ -std=c++20 -c -g -DUSE_DECODER -DUSE_VNS_PURE ../../src/basic.cc -I../../include/ -I${JSON_PATH} -o basic.o

clean:
	rm -r ./bin
	rm *.o
