#include "basic.h"
#include <cstdio>
#include <fstream>

using namespace std;

inline double elapsed_time(hrc::time_point from) {
  return duration_cast<duration<double>>(hrc::now() - from).count();
}

vector<Solution> vns() {
  const int MAX_TIME = 300;
  auto start = high_resolution_clock::now();

  Solution start_sol = CH_VRBSP();
  Solution inc = start_sol;
  printf("%.4lf,%.4lf\n", inc.throughput_, elapsed_time(start));

  Solution inc_20 = convertTo20MhzSol(inc);

  int K_MUL = max(1, N / 100);
  int K_MAX = N;
  while (!stop(start, MAX_TIME)) {
    int k = 1;
    Solution candidate = inc_20;
    while (!stop(start, MAX_TIME) && k <= K_MAX) {
      candidate = perturbation(candidate, k * K_MUL);
      candidate = local_search(candidate);

      if (candidate > inc_20) {
        k = 1;
        inc_20 = candidate;
        printf("%.4lf,%.4lf\n", inc_20.throughput_, elapsed_time(start));
      } else
        k += 1;
    }
  }

#if defined(USE_VNS_PURE)
  Solution best_multiple = multipleRepresentation(inc_20);
  double of = optimal_partitioning_global(best_multiple);

  inc = best_multiple.get_optimal_solution();
  assert(approximatelyEqual(of, inc_20.throughput_));
  auto math_sol = gurobi_sol(inc);

#else
  assert(false);
#endif

  return {start_sol, inc};
}

int main(int argc, char **argv) {
  freopen(argv[1], "r", stdin);
  read_data();

  vector<Solution> ans = vns();

  for (int i = 0; const auto path : {argv[2], argv[3]}) {
    std::ofstream o(path);
    json seila = json(ans[i++]);
    o << seila.dump(4);
    o.close();
  }

  return 0;
}
