#include "DifferentialEvolution.h"
#include "HeuristicDecoder.h"
#include <chrono>

using namespace std;

const double LOWER_BOUND = 0.0;
const double UPPER_BOUND = 1.0;

const double f__ = 0.5;
const double cr__ = 0.5;
int numberOfDifferenceVectors = 1;

int populationSize;
int numberVariables;
int evaluations, maxEvaluations;
int solutionsToSelect;
Solution decoder;
enum STOPPING stopCriteria;
enum CROSSOVER_TYPE crossoverType;
enum MUTATION_TYPE mutationType;
Individual bestSolution;
double maximumTime;
std::chrono::steady_clock::time_point startTime;

enum DE_VARIANT defineVariant(string __variant) {
    if (__variant == "DE_BEST_1_BIN") {
        return BEST_1_BIN;
    } else if (__variant == "DE_BEST_1_EXP") {
        return BEST_1_EXP;
    } else if (__variant == "DE_BEST_2_BIN") {
        return BEST_2_BIN;
    } else if (__variant == "DE_BEST_2_EXP") {
        return BEST_2_EXP;
    } else if (__variant == "DE_RAND_1_BIN") {
        return RAND_1_BIN;
    } else if (__variant == "DE_RAND_1_EXP") {
        return RAND_1_EXP;
    } else if (__variant == "DE_RAND_2_BIN") {
        return RAND_2_BIN;
    } else if (__variant == "DE_RAND_2_EXP") {
        return RAND_2_EXP;
    } else {
        fprintf(stderr, "INVALID VARIANT\n");
        exit(-1);
    }
}

inline void printPopulation(const vector<Individual> aux) {
    for (int i = 0; i < aux.size(); i++) {

        vector<double> variables = aux[i].getVariables();
        double obj = decoder.decode(variables);

        fprintf(stderr, "%d:", i);
        for (int i = 0; i < variables.size(); i++) {
            fprintf(stderr, " %.4lf ", variables[i]);
        }
        fprintf(stderr, "  ==> %.5lf\n", obj);
    }
}

bool isStoppingCriteriaReached() {
    if (stopCriteria == TIMELIMIT) {
        return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - startTime).count() >= maximumTime;
    } else if (stopCriteria == ITERATION) {
        return evaluations > maxEvaluations;
    }

    return false; // TODO
}

vector<Individual> createInitialPopulation() {
    vector<Individual> initialPop;

    for (int i = 0; i < populationSize; i++) {
        vector<double> indVariables;
        for (int i = 0; i < numberVariables; i++) {
            indVariables.emplace_back(rng.rand());
        }

        initialPop.emplace_back(0.0, indVariables);
    }

    return initialPop;
}

void evaluatePopulation(vector<Individual> &population) {
    double MINOBJ = 1000000007;
    for (int i = 0; i < populationSize; i++) {
        vector<double> variables(population[i].getVariables());

        double obj = decoder.decode(variables);

        population[i].setObjective(obj);
        if (obj < MINOBJ) {
            MINOBJ = obj;
            bestSolution = population[i];
        }
    }
}

void initProgress() {
    if (stopCriteria == TIMELIMIT) {
        startTime = std::chrono::high_resolution_clock::now();
    } else if (stopCriteria == ITERATION) {
        evaluations = populationSize;
    }
}

vector<Individual> selection(const vector<Individual> &solutionList) {
    unordered_set<int> toReturn;

    assert(solutionList.size() > 0);

    do {
        int randomIdx = rng.randInt(solutionList.size() - 1);

        while (toReturn.count(randomIdx))
            randomIdx = rng.randInt(solutionList.size() - 1);

        toReturn.insert(randomIdx);
    } while (toReturn.size() < solutionsToSelect);

    vector<Individual> inds;
    for (const int el : toReturn)
        inds.emplace_back(solutionList[el]);

    return inds;
}

double mutate(vector<vector<double>> parent, int index) {
    double value = -1000000007;

    if (mutationType == RAND) {
        if (numberOfDifferenceVectors == 1) {
            return parent[2][index] + f__ * (parent[0][index] - parent[1][index]);
        } else if (numberOfDifferenceVectors == 2) {
            return parent[4][index] + f__ * (parent[0][index] - parent[1][index]) +
                   f__ * (parent[2][index] - parent[3][index]);
        }
    } else if (mutationType == BEST) {
        if (numberOfDifferenceVectors == 1) {
            return bestSolution.getVariableValue(index) +
                   f__ * (parent[0][index] - parent[1][index]);
        } else if (numberOfDifferenceVectors == 2) {
            return bestSolution.getVariableValue(index) +
                   f__ * (parent[0][index] - parent[1][index]) +
                   f__ * (parent[2][index] - parent[3][index]);
        }
    } else if (mutationType == RAND_TO_BEST) {
        // TODO
    }

    return value;
}

vector<Individual> crossover(const vector<Individual> &parentSolutions,
                             Individual child) { // TODO: Maybe there is a bug here
    vector<vector<double>> parent;

    double requiredParents = 1 + numberOfDifferenceVectors * 2;
    assert(requiredParents == solutionsToSelect);

    for (int i = 0; i < requiredParents; i++)
        parent.emplace_back(parentSolutions[i].getVariables());

    int jrand = rng.randInt(numberVariables - 1);

    if (crossoverType == BINARY) {
        for (int j = 0; j < numberVariables; j++) {
            if (rng.rand() || j == jrand) {
                double value = mutate(parent, j);
                child.setVariableValue(j, value);
            }
        }
    } else if (crossoverType == EXPONENTIAL) {
        int j = rng.randInt(numberVariables - 1);
        int l = 0;

        do {
            double value = mutate(parent, j);

            child.setVariableValue(j, value);
            j = (j + 1) % numberVariables;
            l++;
        } while (rng.rand() < cr__ && (l < numberVariables));
    }

    // REPAIR BOUNDS
    vector<double> childVariables = child.getVariables();
    for (int i = 0; i < childVariables.size(); i++) {
        double varValue = childVariables[i];
        if (varValue < LOWER_BOUND) {
            varValue = LOWER_BOUND;
        }

        if (varValue > UPPER_BOUND) {
            varValue = UPPER_BOUND;
        }

        childVariables[i] = varValue;
    }

    child.setVariables(childVariables);
    return {child};
}

vector<Individual> reproduction(vector<Individual> matingPopulation) {
    vector<Individual> offspringPopulation;

    for (int i = 0; i < populationSize; i++) {
        vector<Individual> parents = selection(matingPopulation);
        vector<Individual> children = crossover(parents, matingPopulation[i]);

        offspringPopulation.emplace_back(children[0]);
    }

    return offspringPopulation;
}

vector<Individual> replacement(vector<Individual> &population,
                               vector<Individual> &offspringPopulation) {
    assert(population.size() == offspringPopulation.size());
    vector<Individual> pop;

    for (int i = 0; i < populationSize; i++) {
        if (population[i] < offspringPopulation[i]) {
            pop.emplace_back(population[i]);
        } else {
            pop.emplace_back(offspringPopulation[i]);
        }
    }

    sort(pop.begin(), pop.end());
    bestSolution = pop[0];
    return pop;
}

inline void updateProgress() { evaluations += populationSize; }

vector<Individual> run() {
    vector<Individual> offspringPopulation;

    vector<Individual> _population = createInitialPopulation();
    evaluatePopulation(_population);
    initProgress();

    while (!isStoppingCriteriaReached()) {
        offspringPopulation = reproduction(_population);
        evaluatePopulation(offspringPopulation);
        replacement(_population, offspringPopulation);
        updateProgress();
    }

    sort(_population.begin(), _population.end());
    return _population;
}

// variante, tempo limite, criterio de parada, objectiveFile, solutionFile, instanciaFile
void init(int argc, char **argv, FILE **solutionFile, FILE **objectivesFile) {
    enum DE_VARIANT variant = BEST_1_BIN;

    string path_input = "../instances/md-vrbsp/U_";
    path_input += string(argv[1]) + "/U_";
    // path_input += "/MD-VRBSP_U_";
    path_input += string(argv[1]);
    path_input += "_";
    path_input += string(argv[2]);
    path_input += ".txt";
    if (!path_input.empty()) {
        fprintf(stderr, "trying to open input file %s\n", path_input.c_str());
        freopen(path_input.c_str(), "r", stdin);
    }
    
    maximumTime = stoi(argv[4]);

    if (stdin == nullptr) {
        fprintf(stderr, "error opening input file (stdin)\n");
        exit(13);
    }

    loadData();
    populationSize = 100;
    numberVariables = 2 * nConnections;
    // numberVariables = nConnections;    

    string solFile = string(string(argv[3])  + string("/solutionFile_de") + string(argv[2]) + string(".txt"));      
    *solutionFile = fopen(solFile.c_str(), "w");

    string objFile = string(string(argv[3]) + string("/objectives_de.txt"));
    *objectivesFile = fopen(objFile.c_str(), "a");

    stopCriteria = TIMELIMIT;
    maximumTime = stoi(argv[4]);
    // enum DE_VARIANT variant = defineVariant(argv[5]);

    switch (variant) {
    case RAND_1_BIN:
    case RAND_1_EXP:
    case BEST_1_BIN:
    case BEST_1_EXP:
    case RAND_TO_BEST_1_BIN:
    case RAND_TO_BEST_1_EXP:
    case CURRENT_TO_RAND_1_BIN:
    case CURRENT_TO_RAND_1_EXP:
        numberOfDifferenceVectors = 1;
        break;
    case RAND_2_BIN:
    case RAND_2_EXP:
    case BEST_2_BIN:
    case BEST_2_EXP:
        numberOfDifferenceVectors = 2;
        break;
    default:
        fprintf(stderr, "ERROR IN VARIANT TYPE\n");
        exit(-1);
    }

    switch (variant) {
    case RAND_1_BIN:
    case RAND_1_EXP:
    case BEST_1_BIN:
    case BEST_1_EXP:
    case RAND_TO_BEST_1_BIN:
    case RAND_TO_BEST_1_EXP:
    case CURRENT_TO_RAND_1_BIN:
    case CURRENT_TO_RAND_1_EXP:
        solutionsToSelect = 3;
        break;
    case RAND_2_BIN:
    case RAND_2_EXP:
    case BEST_2_BIN:
    case BEST_2_EXP:
        solutionsToSelect = 5;
        break;
    default:
        fprintf(stderr, "ERROR IN solutionsToSelect\n");
        exit(-1);
    }

    switch (variant) {
    case RAND_1_BIN:
    case BEST_1_BIN:
    case RAND_TO_BEST_1_BIN:
    case CURRENT_TO_RAND_1_BIN:
    case RAND_2_BIN:
    case BEST_2_BIN:
        crossoverType = BINARY;
        break;
    case RAND_1_EXP:
    case BEST_1_EXP:
    case RAND_TO_BEST_1_EXP:
    case CURRENT_TO_RAND_1_EXP:
    case RAND_2_EXP:
    case BEST_2_EXP:
        crossoverType = EXPONENTIAL;
        break;
    default:
        fprintf(stderr, "ERROR IN CROSSOVER VARIANT\n");
        exit(-1);
    }

    switch (variant) {
    case RAND_1_BIN:
    case RAND_1_EXP:
    case RAND_2_BIN:
    case RAND_2_EXP:
        mutationType = RAND;
        break;
    case BEST_1_BIN:
    case BEST_1_EXP:
    case BEST_2_BIN:
    case BEST_2_EXP:
        mutationType = BEST;
        break;
    case CURRENT_TO_RAND_1_BIN:
    case CURRENT_TO_RAND_1_EXP:
        mutationType = CURRENT_TO_RAND;
        break;
    case RAND_TO_BEST_1_BIN:
    case RAND_TO_BEST_1_EXP:
        mutationType = RAND_TO_BEST;
        break;
    default:
        fprintf(stderr, "ERROR IN MUTATION VARIANT\n");
        exit(-1);
    }

    if (stdin == nullptr) {
        fprintf(stderr, "error opening input file (stdin)\n");
        exit(13);
    }

    fprintf(stdout, "[DIFFERENTIAL EVOLUTION variant %s] will execute for %lf seconds\n", argv[5],
            maximumTime);
}

int main(int argc, char *argv[]) { // DE_RAND_1_BIN
    FILE *solutionFile = nullptr, *objectivesFile = nullptr;
    init(argc, argv, &solutionFile, &objectivesFile);

    vector<Individual> population = run();
    sort(population.begin(), population.end());

    if (solutionFile != nullptr) {
        vector<double> ans = bestSolution.getVariables();
        for (int i = 0; i < ans.size(); i++) {
            fprintf(solutionFile, "%lf ", ans[i]);
        }
        fprintf(solutionFile, "\n");
    } else {
        fprintf(stderr, "solutionFile is null!\n");
        exit(13);
    }

    if (objectivesFile != nullptr) {
        fprintf(objectivesFile, "%lf %d\n", bestSolution.getObjective(), evaluations);
    } else {
        fprintf(stderr, "objectivesFiles is null!\n");
        exit(13);
    }

    fclose(solutionFile);
    fclose(objectivesFile);
    return 0;
}
