brkga: decoder.o
	g++ -std=c++11 -c -g HeuristicBRKGA.cpp
	g++ -std=c++11 -o decoder_brkga HeuristicBRKGA.o HeuristicDecoder.o
#	rm *.o

de: decoder.o
	g++ -std=c++11 -c -g DifferentialEvolution.cpp
	g++ -std=c++11 -o decoder_de DifferentialEvolution.o HeuristicDecoder.o

decoder.o:
	g++ -std=c++11 -c -g HeuristicDecoder.cpp

