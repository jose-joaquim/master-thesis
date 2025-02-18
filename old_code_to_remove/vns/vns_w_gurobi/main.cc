// run 128 1 results/vns-vrbsp/U_128 10
#include "common.h"

using namespace std;

Solution vns(string filePrefix) {
  // FILE *outFile = nullptr;
  // if (!filePrefix.empty()) outFile = fopen(filePrefix.c_str(), "w");

  auto ch_start = high_resolution_clock::now();
  Solution inc = CH_VRBSP();

  inc(0).spectrums.emplace_back(Spectrum(0, 0, {Channel(0)}));

  int K_MUL = max(1, N / 100);
  int K_MAX = 10;
  double currTime = 0.0;
  start = high_resolution_clock::now();
  while (!stop()) {
    int k = 1;
    while (!stop() && k <= K_MAX) {
      Solution candidate = perturbation(inc, k * K_MUL);
      candidate = ls_math_optimizer(candidate);

      if (compareObjectives(candidate, inc) < 0) {
        k = 1;
        inc = candidate;
      } else
        k += 1;
    }
  }

  return inc;
}

int main(int argc, char **argv) {
  freopen(argv[3], "r", stdin);
  read_data();

  maximumTime = 300;
  string objective_improvements = string(argv[3]);
  objective_improvements += "/objective_improvements" + string(argv[2]);
  objective_improvements += ".txt";

  Solution aux = vns(objective_improvements);

  // string path_out = string(argv[3]);
  // path_out += "/solution" + string(argv[2]);
  // path_out += ".txt";
  // cout << path_out << endl;
  // FILE *file_out = fopen(path_out.c_str(), "w");
  // print_solution_to_file(aux, file_out);
  return 0;
}
