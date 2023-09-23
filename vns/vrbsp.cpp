// run 128 1 results/vns-vrbsp/U_128 10
#include "common.h"

using namespace std;

Solution vns(string filePrefix) {
  // FILE *outFile = nullptr;
  // if (!filePrefix.empty()) outFile = fopen(filePrefix.c_str(), "w");

  auto ch_start = high_resolution_clock::now();
  Solution incumbent = CH_VRBSP();

  // auto elapsed =
  //     duration_cast<duration<double>>(high_resolution_clock::now() -
  //     ch_start)
  //         .count();
  // fprintf(outFile, "%.3lf %.3lf\n", elapsed, incumbent.throughput);
  // printf("%.3lf %.3lf\n", elapsed, incumbent.throughput);

  Solution delta = convertTo20MhzSol(incumbent);
  Solution local_max = delta;

  int K_MUL = max(1, n_connections / 100);
  int K_MAX = 10;
  double currTime = 0.0;
  start = high_resolution_clock::now();
  while (!stop()) {
    int k = 1;
    while (!stop() && k <= K_MAX) {
      delta = local_max;
      perturbation(delta, k * K_MUL);

      Solution multiple = multipleRepresentation(delta);
      setDP(multiple);
      delta.throughput = calcDP(multiple);

      exit(10);

      Solution explicit_sol = local_search(multiple, delta);
      fix_channels(explicit_sol);

      delta = convertTo20MhzSol(explicit_sol);

      if (compareObjectives(delta, local_max) < 0) {
        k = 1;
        local_max = delta;
      } else
        k += 1;

      if (compareObjectives(local_max, incumbent) < 0) {
        currTime = duration_cast<duration<double>>(
                       high_resolution_clock::now() - start)
                       .count();
        // if (outFile != nullptr)
        //   fprintf(outFile, "%.3lf %.3lf\n", currTime, local_max.throughput);

        printf("%.3lf %.3lf %.3lf\n", currTime, incumbent.throughput,
               local_max.throughput);
        incumbent = explicit_sol;
      }
    }
  }

  return incumbent;
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
