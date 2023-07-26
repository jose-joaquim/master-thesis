#include "mdvrbsp.h"

#include "common.h"

double computeViolation(const Solution &sol) {
  //   if (of_opt == MINMAX) {
  //     sol.violation = 0.0;
  //     for (TimeSlot &ts : sol.slots)
  //       for (Spectrum &sp : ts.spectrums)
  //         for (Channel &ch : sp.channels) {
  //           double higher = 0.0;
  //           for (Connection &conn : ch.connections) {
  //             conn.violation = max(0.0, GMA[conn.id] - conn.throughput);
  //             higher = max(higher, conn.violation);
  //           }
  //
  //           ch.violation = higher;
  //           sol.violation = max(sol.violation, ch.violation);
  //         }
  //   } else {
  //     sol.violation = 0.0;
  //     for (TimeSlot &ts : sol.slots)
  //       for (Spectrum &sp : ts.spectrums)
  //         for (Channel &ch : sp.channels) {
  //           double sum = 0.0;
  //           for (Connection &conn : ch.connections) {
  //             conn.violation = max(0.0, GMA[conn.id] - conn.throughput);
  //             sum += max(0.0, conn.violation);
  //           }
  //
  //           ch.violation = sum;
  //           sol.violation += ch.violation;
  //         }
  //   }
  //
  //   return sol.violation;
  return 0.0;  // TODO
}

bool yFeasible(Solution &sol) {
  double violation = computeViolation(sol);
  return essentiallyEqual(violation, 0.0);
}

Solution DTS(const Solution sol) { return {}; }

Solution VNS(Solution sol) { return {}; }

bool mdvrbsp_feasible(const Solution &ret) {
  int cnt = 0;
  set<int> scheduled;
  for (const TimeSlot &ts : ret.slots) {
    for (const Spectrum &sp : ts.spectrums) {
      for (const Channel &ch : sp.channels) {
        for (const Connection &conn : ch.connections) {
          // printf("conn %d\n", conn.id);
          cnt += 1;
          scheduled.insert(conn.id);
          if (definitelyLessThan(conn.throughput, GMA[conn.id])) {
            printf("ops... conn %d %lf %lf\n", conn.id, conn.throughput,
                   GMA[conn.id]);
            return false;
          }
        }
      }
    }
  }

  if (cnt != n_connections) {
    for (auto x : scheduled) printf("%d ", x);

    puts("");
  }

  return cnt >= n_connections;
}

Solution reductionHeuristic(char **argv) {
  // Solution S_star = CH();
  //
  // if (S_star.slots.size() == 1) {
  //   assert(yFeasible(S_star));
  //   printf(
  //       "Constructive heuristic found feasible solution with only one "
  //       "time-slots.\n");
  //   return S_star;
  // }
  //
  // printf("initial solution has %lu time-slots\n", S_star.slots.size());
  // while (!stop() && S_star.slots.size() > 1) {
  //   Solution S1 = DTS(S_star);
  //
  //   int cnt = 1;
  //   while (yFeasible(S1) && ++cnt && S1.slots.size() > 1LU) {
  //     S_star = S1;
  //     S1 = DTS(S1);
  //   }
  //
  //   if (yFeasible(S1) && S1.slots.size() == 1LU) return S1;
  //
  //   printf("removed %d ts. Violation is now %lf\n", cnt, S1.violation);
  //   S1 = VNS(S1);
  //
  //   if (yFeasible(S1)) {
  //     printf("found sol w/ violation %lf and %lu ts\n", S1.violation,
  //            S1.slots.size());
  //     S_star = S1;
  //   } else
  //     puts("vns did not found feasible solution");
  // }
  //
  // return S_star;

  return Solution();
}

optional<Channel> try_insert(int conn_id, Channel ch) {
  ch = insertInChannel(ch, conn_id);

  for (Connection &conn : ch.connections)
    if (definitelyLessThan(conn.throughput, GMA[conn.id])) return nullopt;

  return ch;
}

bool can_split(const Channel &ch) {
  if (ch.bandwidth <= 20) return false;

  for (const auto &conn : ch.connections)
    if (definitelyLessThan(dataRates[11][bToIdx(ch.bandwidth)], GMA[conn.id]))
      return false;

  Channel aux(ch.bandwidth / 2);
  for (auto &conn : ch.connections)
    if (try_insert(conn.id, aux) == nullopt) return false;

  return true;
}

Solution CH() {
  TimeSlot dummy_ts(init_conf);
  Solution ret({dummy_ts});

  vector<int> links(n_connections);
  iota(links.begin(), links.end(), 0);
  shuffle(links.begin(), links.end(), whatever);

  for (int conn : links) {
    for (int t = 0; t < ret.slots.size(); t++) {
      for (int s = 0; s < ret(t).spectrums.size(); s++) {
        for (int c = 0; c < ret(t, s).channels.size(); c++) {
          if (make_tuple(t, s, c) == zeroChannel) continue;

          auto rst = try_insert(conn, ret(t, s, c));

          if (rst) {
            ret(t, s, c) = *rst;
            goto NEXT_CONN;
          }

          if (can_split(ret(t, s, c))) {
            vector<Channel> ch_split = split(ret(t, s, c));

            for (int i = 0; i < ch_split.size(); ++i) {
              auto cp_channel = try_insert(conn, ch_split[i]);

              if (cp_channel) {
                swap(ret(t, s, c), ret(t, s).channels.back());
                ret(t, s).channels.pop_back();
                ret(t, s).channels.emplace_back(*cp_channel);
                ret(t, s).channels.emplace_back(ch_split[i == 0]);
                goto NEXT_CONN;
              }
            }
          }
        }
      }
    }

    {
      TimeSlot new_ts{dummy_ts};
      for (int s = 0; s < new_ts.spectrums.size(); s++) {
        for (int c = 0; c < new_ts(s).channels.size(); c++) {
          int bw = new_ts(s, c).bandwidth;
          if (bw >= 160) {
            Connection new_conn{conn};
            new_conn.SINR = 1000000007;
            new_conn.throughput = dataRates[11][bwIdx(bw)];
            new_ts(s, c).connections.emplace_back(new_conn);
            ret.slots.emplace_back(new_ts);
            goto NEXT_CONN;
          }
        }
      }
    }

  NEXT_CONN:;
    // computeThroughput(ret);
    // assert(mdvrbsp_feasible(ret));
  }

  vector<Channel> aux;
  aux.emplace_back(Channel());
  ret.slots[0].spectrums.emplace_back(0, 0, aux);

  computeThroughput(ret);

  assert(mdvrbsp_feasible(ret));
  return ret;
}

int main(int argc, char **argv) {
  freopen(argv[3], "r", stdin);
  read_data(GMA);

  Solution inc = CH();
  cout << inc.slots.size() << endl;

  string fsol_s = "sol" + string(argv[1]) + "_" + string(argv[2]) + ".sol";
  FILE *fileSol = fopen(fsol_s.c_str(), "w");
  print_solution_to_gurobi(inc, fileSol);
  fclose(fileSol);
  return 0;
}
