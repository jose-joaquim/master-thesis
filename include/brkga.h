/*
 * BRKGA.h
 *
 * This class encapsulates a Biased Random-key Genetic Algorithm (for
 * minimization problems) with K independent Populations stored in two vectors
 * of Population, current and previous. It supports multi-threading via OpenMP,
 * and implements the following key methods:
 *
 * - BRKGA() constructor: initializes the populations with parameters described
 * below.
 * - evolve() operator: evolve each Population following the BRKGA methodology.
 * This method supports OpenMP to evolve up to K independent Populations in
 * parallel. Please note that double Decoder::decode(...) MUST be thread-safe.
 *
 * Required hyperparameters:
 * - n: number of genes in each chromosome
 * - p: number of elements in each population
 * - pe: pct of elite items into each population
 * - pm: pct of mutants introduced at each generation into the population
 * - rhoe: probability that an offspring inherits the allele of its elite parent
 *
 * Optional parameters:
 * - K: number of independent Populations
 * - MAX_THREADS: number of threads to perform parallel decoding -- WARNING:
 * Decoder::decode() MUST be thread-safe!
 *
 * Required templates are:
 * RNG: random number generator that implements the methods below.
 *     - RNG(unsigned long seed) to initialize a new RNG with 'seed'
 *     - double rand() to return a double precision random deviate in range
 * [0,1)
 *     - unsigned long randInt() to return a >=32-bit unsigned random deviate in
 * range [0,2^32-1)
 *     - unsigned long randInt(N) to return a unsigned random deviate in range
 * [0, N] with N < 2^32
 *
 * Decoder: problem-specific decoder that implements any of the decode methods
 * outlined below. When compiling and linking BRKGA with -fopenmp (i.e., with
 * multithreading support via OpenMP), the method must be thread-safe.
 *     - double decode(const vector< double >& chromosome) const, if you don't
 * want to change chromosomes inside the framework, or
 *     - double decode(vector< double >& chromosome) const, if you'd like to
 * update a chromosome
 *
 *  Created on : Jun 22, 2010 by rtoso
 *  Last update: Sep 28, 2010 by rtoso
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 */

#ifndef BRKGA_H
#define BRKGA_H

#include "MTRand.h"
#include "population.h"
#include <algorithm>
#include <exception>
#include <omp.h>
#include <stdexcept>

extern int evaluations;

template <class Decoder, class RNG> class BRKGA {
public:
  /*
   * Default constructor
   * Required hyperparameters:
   * - n: number of genes in each chromosome
   * - p: number of elements in each population
   * - pe: pct of elite items into each population
   * - pm: pct of mutants introduced at each generation into the population
   * - rhoe: probability that an offspring inherits the allele of its elite
   * parent
   *
   * Optional parameters:
   * - K: number of independent Populations
   * - MAX_THREADS: number of threads to perform parallel decoding
   *                WARNING: Decoder::decode() MUST be thread-safe; safe if
   * implemented as
   *                + double Decoder::decode(std::vector< double >& chromosome)
   * const
   */
  BRKGA(unsigned n, unsigned p, double pe, double pm, double rhoe,
        const Decoder &refDecoder, RNG &refRNG, unsigned K = 1,
        unsigned MAX_THREADS = 1);

  /**
   * Destructor
   */
  ~BRKGA();

  /**
   * Resets all populations with brand new keys
   */
  void reset();

  /**
   * Evolve the current populations following the guidelines of BRKGAs
   * @param generations number of generations (must be even and nonzero)
   * @param J interval to exchange elite chromosomes (must be even; 0 ==> no
   * synchronization)
   * @param M number of elite chromosomes to select from each population in
   * order to exchange
   */
  void evolve(unsigned generations = 1);

  /**
   * Exchange elite-solutions between the populations
   * @param M number of elite chromosomes to select from each population
   */
  void exchangeElite(unsigned M);

  /**
   * Returns the current population
   */
  const Population &getPopulation(unsigned k = 0) const;

  /**
   * Returns the chromosome with best fitness so far among all populations
   */
  const std::vector<double> &getBestChromosome() const;

  /**
   * Returns the best fitness found so far among all populations
   */
  double getBestFitness() const;

  // Return copies to the internal parameters:
  unsigned getN() const;
  unsigned getP() const;
  unsigned getPe() const;
  unsigned getPm() const;
  unsigned getPo() const;
  double getRhoe() const;
  unsigned getK() const;
  unsigned getMAX_THREADS() const;

private:
  // Hyperparameters:
  const unsigned n;  // number of genes in the chromosome
  const unsigned p;  // number of elements in the population
  const unsigned pe; // number of elite items in the population
  const unsigned
      pm; // number of mutants introduced at each generation into the population
  const double rhoe; // probability that an offspring inherits the allele of its
                     // elite parent

  // Templates:
  RNG &refRNG;               // reference to the random number generator
  const Decoder &refDecoder; // reference to the problem-dependent Decoder

  // Parallel populations parameters:
  const unsigned K;           // number of independent parallel populations
  const unsigned MAX_THREADS; // number of threads for parallel decoding

  // Data:
  std::vector<Population *> previous; // previous populations
  std::vector<Population *> current;  // current populations

  // Local operations:
  void initialize(
      const unsigned i); // initialize current population 'i' with random keys
  void evolution(Population &curr, Population &next);
  bool isRepeated(const std::vector<double> &chrA,
                  const std::vector<double> &chrB) const;
};

#endif
