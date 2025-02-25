#include "basic.h"
#include <cstdio>

using namespace std;

inline double elapsed_time(hrc::time_point from) {
  return duration_cast<duration<double>>(hrc::now() - from).count();
}

Solution vns() {
  progress_file = fopen("progress.csv", "w");
  auto start = high_resolution_clock::now();
  Solution inc = CH_VRBSP();
  fprintf(progress_file, "%.4lf,%.4lf\n", inc.throughput_, elapsed_time(start));
  Solution inc_20 = convertTo20MhzSol(inc);

  int K_MUL = max(1, N / 100);
  int K_MAX = N;
  double currTime = 0.0;
  while (!stop(start, 30)) {
    int k = 1;
    Solution candidate = inc_20;
    // assert(!inc_20.get_scheduled_connections().empty());
    while (!stop(start, 30) && k <= K_MAX) {
      // assert(!candidate.get_scheduled_connections().empty());
      candidate = perturbation(candidate, k * K_MUL);
      candidate = local_search(candidate);

      if (candidate > inc_20) {
        k = 1;
        inc_20 = candidate;
        // Solution best_multiple = multipleRepresentation(inc_20);
        // double of = optimal_partitioning_global(best_multiple);
        // assert(approximatelyEqual(of, inc_20.throughput_));

        printf("%.4lf,%.4lf\n", inc_20.throughput_, elapsed_time(start));
        fprintf(progress_file, "%.4lf,%.4lf\n", inc_20.throughput_,
                elapsed_time(start));
      } else
        k += 1;
    }
  }

  fclose(progress_file);
  Solution best_multiple = multipleRepresentation(inc_20);
  double of = optimal_partitioning_global(best_multiple);

  assert(approximatelyEqual(of, inc_20.throughput_));
  inc = best_multiple.get_optimal_solution();
  auto math_sol = gurobi_sol(inc);

  return inc;
}

int main(int argc, char **argv) {
  freopen(argv[1], "r", stdin);
  read_data();

  Solution ans = vns();
  return 0;
}
