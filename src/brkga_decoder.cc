#include "MTRand.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/globals.h"
#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "basic.h"
#include "brkga.h"
#include <fstream>

namespace fs = std::filesystem;

int numberVariables;

ABSL_FLAG(
    int, type, 0,
    "Decoder/heuristic type: Classic - 0 | Classic DP - 1 | Fixed DP - 2");
ABSL_FLAG(double, maximum_time, 0.0, "Maximum execution time (seconds).");

ABSL_FLAG(bool, first, true,
          "Whether this is the first run (controls one-time output "
          "directory cleanup).");

ABSL_FLAG(unsigned, p, 100, "Population size.");
ABSL_FLAG(double, pe, 0.05, "Elite population fraction.");
ABSL_FLAG(double, pm, 0.25, "Mutant (bottom) population fraction.");
ABSL_FLAG(double, rhoe, 0.70,
          "Probability of inheriting an allele from the elite parent.");
ABSL_FLAG(unsigned, k, 1, "Number of independent populations.");
ABSL_FLAG(unsigned, maxt, static_cast<unsigned>(omp_get_max_threads()),
          "Maximum number of threads.");
ABSL_FLAG(unsigned, x_intvl, 100, "Exchange interval (generations).");
ABSL_FLAG(unsigned, x_number, 2, "Number of individuals exchanged.");
ABSL_FLAG(unsigned, max_gens, 1000, "Maximum number of generations.");

// type 6: se true, aplica localSearchType6 nas pe elites a cada geracao
// (alem da LS final na melhor solucao). Edite e recompile para alternar.
ABSL_FLAG(bool, use_elite_ls, true,
          "Type 6 only: apply localSearchType6 to elite individuals every "
          "generation (in addition to the final LS on the best solution).");

int type;
double maximumTime;

bool first;

unsigned p;
double pe;
double pm;
double rhoe;
unsigned K;
unsigned MAXT;
unsigned X_INTVL;
unsigned X_NUMBER;
unsigned MAX_GENS;

bool USE_ELITE_LS;

void init(const std::string &instancePath, FILE **solutionFile,
          FILE **objectivesFile, FILE **timeFile, FILE **populationFile,
          FILE **channelsFile) {
  LOG(INFO) << "Loading instance from " << instancePath;
  std::ifstream input_json_path(instancePath);
  json input_json;
  input_json_path >> input_json;
  Instance inst;
  from_json(input_json, inst);

  numberVariables = 2 * inst.N;

  string outputDir;

  ostringstream oss;
  oss << (pe * 100);
  string aux = oss.str();
  if (type == 0)
    outputDir = "../output/output_" + aux + '/' + to_string(inst.N);
  else if (type == 1 || type == 4)
    outputDir = "../output/output_dp_" + aux + '/' + to_string(inst.N);
  else if (type == 3)
    outputDir = "../output/output_random_dp_" + aux + '/' + to_string(inst.N);
  else if (type == 5)
    outputDir =
        "../output/output_new_heuristic_2_" + aux + '/' + to_string(inst.N);
  else if (type == 6)
    outputDir =
        "../output/output_local_search_" + aux + '/' + to_string(inst.N);

  try {
    if (fs::exists(outputDir) && first) {
      fs::remove_all(outputDir);
      first = false;
    }
  } catch (const fs::filesystem_error &e) {
    LOG(FATAL) << "Error cleaning up output directories: " << e.what();
  }

  try {
    fs::create_directories(outputDir);
  } catch (const fs::filesystem_error &e) {
    fprintf(stderr, "Error creating directory %s: %s\n", outputDir.c_str(),
            e.what());
    exit(1);
  }

  string solFile = outputDir + "/solution.txt";
  *solutionFile = fopen(solFile.c_str(), "a");

  string objFile = outputDir + "/objectives.txt";
  *objectivesFile = fopen(objFile.c_str(), "a");

  string tFile = outputDir + "/time.txt";
  *timeFile = fopen(tFile.c_str(), "a");

  string popFile = outputDir + "/population.txt";
  *populationFile = fopen(popFile.c_str(), "a");

  // Composicao de canais da solucao final: uma linha por instancia com
  // "c20 c40 c80 c160".
  string chanFile = outputDir + "/channels.txt";
  *channelsFile = fopen(chanFile.c_str(), "a");
}

int main(int argc, char **argv) {
  absl::SetStderrThreshold(absl::LogSeverityAtLeast::kInfo);
  absl::InitializeLog();
  absl::ParseCommandLine(argc, argv);

  type = absl::GetFlag(FLAGS_type);
  maximumTime = absl::GetFlag(FLAGS_maximum_time);

  first = absl::GetFlag(FLAGS_first);

  p = absl::GetFlag(FLAGS_p);
  pe = absl::GetFlag(FLAGS_pe);
  pm = absl::GetFlag(FLAGS_pm);
  rhoe = absl::GetFlag(FLAGS_rhoe);
  K = absl::GetFlag(FLAGS_k);
  MAXT = absl::GetFlag(FLAGS_maxt);
  X_INTVL = absl::GetFlag(FLAGS_x_intvl);
  X_NUMBER = absl::GetFlag(FLAGS_x_number);
  MAX_GENS = absl::GetFlag(FLAGS_max_gens);

  USE_ELITE_LS = absl::GetFlag(FLAGS_use_elite_ls);

  // So faz sentido para o type 6 (unico decoder com localSearchType6).
  bool applyEliteLS = (type == 6) && USE_ELITE_LS;
  FILE *solutionFile = nullptr, *objectivesFile = nullptr, *timeFile = nullptr,
       *populationFile = nullptr, *channelsFile = nullptr;

  // evaluations = 0;
  const unsigned n = numberVariables;
  Solution decoder;
  MTRand rng;
  BRKGA<Solution, MTRand> algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
  double TempoExecTotal = 0.0, TempoFO_Star = 0.0, FO_Star = 1000000007,
         FO_Min = -1000000007;
  int bestGeneration = 0, minGeneration = 0;
  int iterSemMelhora, iterMax = 10, quantIteracoes = 0, bestIteration = 0;

  clock_t TempoFO_StarInic;

  TempoFO_StarInic = clock();

  bestGeneration = 0;
  minGeneration = 0;

  iterSemMelhora = 0;
  iterMax = 1;

  quantIteracoes = 0;
  bestIteration = 0;

  unsigned generation = 0;
  const auto DURATION_LIMIT = std::chrono::minutes(5);
  for (int i = 0; i < 1000; i++) {
    LOG_IF(INFO, i % 50) << "AQUI";

    algorithm.evolve();

    if (populationFile != nullptr) {
      for (unsigned k = 0; k < p; k++) {
        // double val = -1.0 * algorithm.getPopulationFitness(0, k);
        // fprintf(populationFile, "%lf ", val);
      }
      // fprintf(populationFile, "\n"); // Nova linha = Nova geração
    }

    if ((++generation) % X_INTVL == 0) {
      algorithm.exchangeElite(X_NUMBER);
    }

    if (algorithm.getBestFitness() < FO_Star) {
      TempoFO_Star = (((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);
      FO_Star = algorithm.getBestFitness();
      bestGeneration = generation;
      bestIteration = quantIteracoes;
    }
  }

  TempoExecTotal = (((double)(clock() - TempoFO_StarInic)) / CLOCKS_PER_SEC);

  if (solutionFile != nullptr) {
    vector<double> best = algorithm.getBestChromosome();
    for (int i = 0; i < best.size(); i++)
      fprintf(solutionFile, "%lf ", best[i]);
    fprintf(solutionFile, "\n");
  } else {
    cout << stderr << "solutionFile is null!" << endl;
    exit(13);
  }

  // type 6 (Pacote 1): aplica a busca local SÓ na melhor solução, ao final.
  // Aproveita esse decode da melhor solução para coletar a composição de canais
  // (channels.txt). Para os demais tipos, decodifica a melhor uma vez só p/ os
  // canais.
  bool collectChannels, applyFinalLS;
  long long count_improvements = 0;
  map<int, long long> channelBandCount;
  double finalFitness = algorithm.getBestFitness();
  {
    vector<double> bestChromo = algorithm.getBestChromosome();
    collectChannels = true;
    if (type == 6)
      applyFinalLS = true;
    double f = decoder.decode(bestChromo); // preenche channelBandCount
    if (type == 6) {
      finalFitness = f;
      applyFinalLS = false;
    }
    collectChannels = false;
  }

  // if (objectivesFile != nullptr)
  //   fprintf(objectivesFile, "%lf %lld %lld\n", -1.0 * finalFitness,
  //   evaluations,
  //           count_improvements);
  // else {
  //   cout << stderr << "objectivesFiles is null!" << endl;
  //   exit(13);
  // }

  // Uma linha por instancia: contagem de canais de 20/40/80/160 MHz na solucao
  // final.
  if (channelsFile != nullptr) {
    fprintf(channelsFile, "%lld %lld %lld %lld\n", channelBandCount[20],
            channelBandCount[40], channelBandCount[80], channelBandCount[160]);
  } else {
    cout << stderr << "channelsFile is null!" << endl;
    exit(13);
  }

  if (timeFile != nullptr)
    fprintf(timeFile, "%lf\n", TempoExecTotal);
  else {
    cout << stderr << "timeFile is null!" << endl;
    exit(13);
  }

  fclose(solutionFile);
  fclose(objectivesFile);
  fclose(timeFile);
  fclose(populationFile);
  fclose(channelsFile);

  return 0;
}
