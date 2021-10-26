#ifndef BRKGA_FF_BEST_DIFFERENTIALEVOLUTION_H
#define BRKGA_FF_BEST_DIFFERENTIALEVOLUTION_H

#include <vector>
#include <unordered_set>
#include "BRKGA.h"

extern const double LOWER_BOUND;
extern const double UPPER_BOUND;

extern const double f__;
extern const double cr__;
extern int numberOfDifferenceVectors;

class Individual {
    double objective;
    std::vector<double> variables;

public:

    Individual() {
        objective = 0.0;
    }

    Individual(const int _objective, const std::vector<double> &_variables) : objective(_objective),
                                                                              variables(_variables) {}

    Individual(std::vector<double> &_variables) : variables(_variables) {
        Individual();
    }

    double getObjective() const {
        return objective;
    }

    void setObjective(double objective) {
        this->objective = objective;
    }

    const std::vector<double> &getVariables() const {
        return variables;
    }

    void setVariables(const std::vector<double> &variables) {
        this->variables = variables;
    }

    void setVariableValue(const int index, const double value) {
        this->variables[index] = value;
    }

    double getVariableValue(const int index) {
        return this->variables[index];
    }

    bool operator<(const Individual &o1) const {
        return objective < o1.objective;
    }
};

enum STOPPING {
    TIMELIMIT,
    ITERATION
};

enum CROSSOVER_TYPE {
    EXPONENTIAL,
    BINARY
};

enum MUTATION_TYPE {
    RAND,
    BEST,
    CURRENT_TO_RAND,
    RAND_TO_BEST
};

enum DE_VARIANT {
    RAND_1_BIN,
    RAND_1_EXP,
    RAND_2_BIN,
    RAND_2_EXP,
    BEST_1_BIN,
    BEST_1_EXP,
    BEST_2_BIN,
    BEST_2_EXP,
    RAND_TO_BEST_1_BIN,
    RAND_TO_BEST_1_EXP,
    CURRENT_TO_RAND_1_BIN,
    CURRENT_TO_RAND_1_EXP
};

std::vector<Individual> run();

bool isStoppingCriteriaReached();

std::vector<Individual> createInitialPopulation();

void initProgress();

std::vector<Individual> selection(const std::vector<Individual> &solutionList);

double mutate(std::vector<std::vector<double> > parent, int index);

std::vector<Individual> crossover(const std::vector<Individual> &parentSolutions, Individual child);

std::vector<Individual> reproduction(std::vector<Individual> matingPopulation);

std::vector<Individual> replacement(std::vector<Individual> &population, std::vector<Individual> &offspringPopulation);

inline void updateProgress();

void evaluatePopulation(std::vector<Individual> &population);

enum DE_VARIANT defineVariant(std::string __variant);

#endif //BRKGA_FF_BEST_DIFFERENTIALEVOLUTION_H
