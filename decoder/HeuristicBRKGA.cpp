#include "BRKGA.h"
#include "HeuristicDecoder.h"
#include "MTRand.h"

using namespace std;

int populationSize;
int numberVariables;
double maximumTime;

// variante, tempo limite, criterio de parada, objectiveFile, solutionFile, instanciaFile
void init(int argc, char **argv, FILE **solutionFile = nullptr, FILE **objectivesFile = nullptr) {
    if (argc != 5) {
        fprintf(stderr,
                "wrong arguments. Provided %d, Must be: stdin, solutionFile, objectiveFile, "
                "timeLimit\n",
                argc);
        exit(13);
    }

    const string openingFile(argv[1]);
    if (!openingFile.empty()) {
        fprintf(stderr, "trying to open input file %s\n", openingFile.c_str());
        freopen(openingFile.c_str(), "r", stdin);
    }

    *solutionFile = fopen(argv[2], "a");
    if (*solutionFile == nullptr) {
        fprintf(stderr, "error opening solutionFile file\n");
        exit(13);
    }

    *objectivesFile = fopen(argv[3], "a");
    if (*objectivesFile == nullptr) {
        fprintf(stderr, "error opening objectivesFile file\n");
        exit(13);
    }

    maximumTime = stoi(argv[4]);

    if (stdin == nullptr) {
        fprintf(stderr, "error opening input file (stdin)\n");
        exit(13);
    }

    loadData();
    populationSize = 100;
    numberVariables = 2 * nConnections;

    fprintf(stdout, "[BRKGA variant %s] will execute for %lf seconds\n", argv[5], maximumTime);
}

int main(int argc, char **argv) {
    FILE *solutionFile = nullptr, *objectivesFile = nullptr;
    init(argc, argv, &solutionFile, &objectivesFile);

    const unsigned p = populationSize; // size of population
    const double pe = 0.25;            // fraction of population to be the elite-set
    const double pm = 0.05;            // fraction of population to be replaced by mutants
    const double rhoe = 0.70; // probability that offspring inherit an allele from elite parent
    const unsigned K =
        1; // number of independent populations         //***Antes era 3 o valor desse parâmetro
    const unsigned MAXT =
        1; // number of threads for parallel decoding   //***Antes era 3 o valor desse parâmetro

    const unsigned X_INTVL = 100;   // exchange best individuals at every 100 generations
    const unsigned X_NUMBER = 2;    // exchange top 2 best
    const unsigned MAX_GENS = 1000; // run for 1000 gens

    string fileName; // TODO

    const unsigned n = numberVariables;

    Solution decoder;
    MTRand rng;
    BRKGA<Solution, MTRand> algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

    double TempoExecTotal = 0.0, TempoFO_Star = 0.0, FO_Star = 1000000007, FO_Min = -1000000007;
    int bestGeneration = 0, minGeneration = 0;
    int iterSemMelhora, iterMax = 10, quantIteracoes = 0, bestIteration = 0;

    clock_t TempoFO_StarInic; // TempoInicial

    TempoFO_StarInic = clock();

    bestGeneration = 0;
    minGeneration = 0;

    iterSemMelhora = 0;
    iterMax = 1;

    quantIteracoes = 0;
    bestIteration = 0;

    unsigned generation = 0; // current generation
    do {
        algorithm.evolve(); // evolve the population for one generation

        if ((++generation) % X_INTVL == 0) {
            algorithm.exchangeElite(X_NUMBER); // exchange top individuals
        }

        if (algorithm.getBestFitness() < FO_Star) {
            TempoFO_Star = (((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);
            FO_Star = algorithm.getBestFitness();
            bestGeneration = generation;
            bestIteration = quantIteracoes;
        }
        quantIteracoes++;
    } while ((((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC) < maximumTime);

    TempoExecTotal = (((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);

    if (solutionFile != nullptr) {
        vector<double> best = algorithm.getBestChromosome();
        for (int i = 0; i < best.size(); i++) {
            fprintf(solutionFile, "%lf ", best[i]);
        }
        fprintf(solutionFile, "\n");
    } else {
        fprintf(stderr, "solutionFile is null!\n");
        exit(13);
    }

    if (objectivesFile != nullptr) {
        fprintf(objectivesFile, "%lf\n", -1.0 * algorithm.getBestFitness());
    } else {
        fprintf(stderr, "objectivesFiles is null!\n");
        exit(13);
    }

    fclose(solutionFile);
    fclose(objectivesFile);
    return 0;
}
