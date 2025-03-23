OMP_PREFIX="/opt/homebrew/opt/libomp"

brkga: decoder.o
	g++ -std=c++20 -c -g -L$OMP_PREFIX/lib -I$OMP_PREFIX/include HeuristicBRKGA.cpp
	g++ -std=c++20 -o ./bin/decoder_brkga HeuristicBRKGA.o HeuristicDecoder.o

decoder.o:
	mkdir -p bin
	g++ -std=c++20 -c -g -DUSE_DECODER ../../src/basic.cc -I../../include/ -o bin/basic.o

