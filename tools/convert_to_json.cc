#include "basic.h"
#include <fstream>
#include <iostream>

// Pull in globals populated by read_data()
extern int N, T;
extern double alfa, NOI, powerSender;
extern double receivers[MAX_CONN][2], senders[MAX_CONN][2];
extern double AFF[MAX_CONN][MAX_CONN];
extern double DM[MAX_CONN][MAX_CONN];
extern vector<vector<double>> DR, SINR, B;
extern vector<double> GMM, BM;
extern vector<vector<int>> overlap;
extern const vector<set<int>> available_channels; // defined in basic.cc

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s out.json < input.txt\n", argv[0]);
    return 1;
  }

  // Parses the legacy plain-text format from stdin into the globals above.
  read_data();

  Instance inst;
  inst.N = N;
  inst.alfa = alfa;
  inst.NOI = NOI;
  inst.powerSender = powerSender;

  inst.receivers.resize(N);
  inst.senders.resize(N);
  for (int i = 0; i < N; i++) {
    inst.receivers[i] = {receivers[i][0], receivers[i][1]};
    inst.senders[i] = {senders[i][0], senders[i][1]};
  }

  // AFF/DM are stored as vector<pair<double,double>> in Instance's current
  // declaration, which only fits a 1-D layout — this is itself a bug to fix
  // alongside the migration (see note below). For now, flatten row-major.
  inst.AFF.clear();
  inst.DM.clear();
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++) {
      inst.AFF.push_back({AFF[i][j], 0.0});
      inst.DM.push_back({DM[i][j], 0.0});
    }

  inst.DR = DR;
  inst.SINR = SINR;
  inst.B = B;
  inst.GMM = GMM;
  inst.BM = BM;
  inst.overlap = overlap;

  inst.C_b.clear();
  for (const auto &s : available_channels)
    inst.C_b.emplace_back(s.begin(), s.end());

  json j = inst;
  std::ofstream out(argv[1]);
  out << j.dump(2);

  return 0;
}
