#include "BRKGA.h"
#include "MTRand.h"
#include "basic.h"

using namespace std;

int populationSize;
int numberVariables;
double maximumTime;

void init(int argc, char **argv, FILE **solutionFile, FILE **objectivesFile) {
    string path_input = argv[1];
    if (!path_input.empty())
        freopen(path_input.c_str(), "r", stdin);

    maximumTime = 10;

    if (stdin == nullptr) {
        fprintf(stderr, "error opening input file (stdin)\n");
        exit(13);
    }

    read_data();
    populationSize = 100;
    numberVariables = 2 * N;
}

void compute_and_save_best_solution(const vector<double> &vars, const double &obj,
                                    char const *const path) {
    vector<pair<double, int>> ranking;
    for (int i = 0; i < vars.size(); i += 2)
        ranking.emplace_back(vars[i], i);

    sort(ranking.begin(), ranking.end());

    vector<int> permutation;
    for (int i = 0; i < ranking.size(); i++)
        permutation.emplace_back(ranking[i].second);

    buildVRBSPSolution(vars, permutation, path);
}

int main(int argc, char **argv) {
    FILE *solutionFile = nullptr, *objectivesFile = nullptr;
    init(argc, argv, &solutionFile, &objectivesFile);

    const unsigned p = populationSize;
    const double pe = 0.25;
    const double pm = 0.05;
    const double rhoe = 0.70;
    const unsigned K = 1;
    const unsigned MAXT = 1;

    const unsigned X_INTVL = 100;
    const unsigned X_NUMBER = 2;

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

    unsigned generation = 0;
    do {
        algorithm.evolve();

        if ((++generation) % X_INTVL == 0) {
            algorithm.exchangeElite(X_NUMBER);
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
    vector<double> best = algorithm.getBestChromosome();
    printf("%lf\n", algorithm.getBestFitness());

    compute_and_save_best_solution(best, algorithm.getBestFitness(), argv[2]);
    return 0;
}
