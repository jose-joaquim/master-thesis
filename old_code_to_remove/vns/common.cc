#include "common.h"

#include <random>

enum CH_PARAMS { RANDOM, GREEDY, MIN_GAMMA, MINIMUM_AFF };
#ifdef MDVRBSP
enum OF_TYPE { MINMAX, SUMVIO };
#endif

int N, n_spectrums;
double alfa, noise, powerSender;
double receivers[MAX_CONN][2], senders[MAX_CONN][2];
double AFF[MAX_CONN][MAX_CONN];
double DM[MAX_CONN][MAX_CONN];
vector<vector<double>> DR, SINR;
vector<int> spectrum_size;

const int MAX_SPECTRUM = 5;
const int MAX_CHANNELS = 45;
const int MAX_SLOTS = 2048;
//  double startTime;
double maximumTime;

high_resolution_clock::time_point start;

int parent[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
int child[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS][2];
double chanOF[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
bool inSolution[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
CH_PARAMS ch_opt;

#ifdef MDVRBSP
OF_TYPE of_opt;
vector<double> gma;
#endif

random_device rd;
mt19937 gen{rd()};
uniform_int_distribution<> dist{0, MAX_CONN};
ti3 zeroChannel = {0, 3, 0};

int overlap[45][45] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0},
    {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
     1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
     0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1}};

bool approximatelyEqual(double a, double b, double epsilon) {
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentiallyEqual(double a, double b, double epsilon) {
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyGreaterThan(double a, double b, double epsilon) {
  return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(double a, double b, double epsilon) {
  return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool stop(void) {
  return duration_cast<duration<double>>(high_resolution_clock::now() - start)
             .count() >= maximumTime;
}

double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
  return hypot((X_si - X_ri), (Y_si - Y_ri));
}

int random_number(int max_num = MAX_CONN) { return dist(gen) % max_num; }

void K_AddDrop(Solution &sol, int K) {
  K = min(K, int(sol(0, 3, 0).connections.size()));
  ti3 channelFrom = {0, 3, 0};
  vector<int> inserted, sizes;
  for (int i = 0; i < K; i++) {
    int idx = random_number(sol(0, 3, 0).connections.size());
    Connection conn = Connection(sol(0, 3, 0, idx));

    int t = random_number(sol.slots.size());
    int a = random_number(sol(t).spectrums.size());
    int b = random_number(sol(t, a).channels.size());
    ti3 channelTo = {t, a, b};

    reinsert(sol, conn, channelFrom, channelTo, true);
  }
}

void K_RemoveAndInserts(Solution &sol, int K) {
  int k = 0;

  while (k < K) {
    int t = random_number(sol.slots.size());
    assert(t >= 0);
    int a = random_number(sol(t).spectrums.size());

    assert(a >= 0);
    int b = random_number(sol(t, a).channels.size());

    assert(b >= 0);

    if (sol(t, a, b).connections.empty()) continue;

    Channel &ch = sol(t, a, b);
    int z = random_number(ch.connections.size());
    Connection conn = ch.connections[z];

    reinsert(sol, conn, make_tuple(t, a, b), zeroChannel, true);
    k++;
  }

  assert(!sol(0, 3, 0).connections.empty());
  K_AddDrop(sol, K);
}

bool allChannels20MHz(const Solution &sol) {
  for (const TimeSlot &ts : sol.slots)
    for (const Spectrum &spectrum : ts.spectrums)
      for (const Channel &channel : spectrum.channels)
        if (channel.bandwidth > 20) return false;

  return true;
}

bool fix_channels(Solution &sol) {
  bool improved = false;
  do {
    improved = false;
    for (int t = 0; t < sol.slots.size(); t++) {
      for (int s = 0; s < sol(t).spectrums.size(); s++) {
        for (int c = 0; c < sol(t, s).channels.size(); c++) {
          if (make_tuple(t, s, c) == zeroChannel) continue;

          int idx = 0;
          while (idx < sol(t, s, c).connections.size()) {
            int bandwidth = sol(t, s, c).bandwidth;
            Connection &conn = sol(t, s, c).connections[idx];
            double aux = computeConnectionThroughput(conn, bandwidth);

            if (approximatelyEqual(aux, 0.0, EPS)) {
#ifdef MDVRBSP
              assert(false);
#endif
              reinsert(sol, conn, make_tuple(t, s, c), zeroChannel, true);
              improved = true;
            } else {
              idx++;
            }
          }
        }
      }
    }
  } while (improved);
  return improved;
}

Solution convertTo20MhzSol(Solution exps) {
  for (TimeSlot &ts : exps.slots)
    for (Spectrum &sp : ts.spectrums)
      sort(sp.channels.rbegin(), sp.channels.rend());

  while (!allChannels20MHz(exps))
    for (TimeSlot &ts : exps.slots)
      for (Spectrum &sp : ts.spectrums)
        for (int c = 0; c < sp.channels.size(); c++) {
          if (sp.channels[c].bandwidth <= 20) continue;

          int newBw = sp(c).bandwidth / 2;
          Channel child1(newBw), child2(newBw);

          vector<Connection> connections = sp(c).connections;

          for (Connection &conn : connections)
            if (dist(gen) % 2)
              child1.connections.emplace_back(conn);
            else
              child2.connections.emplace_back(conn);

          sp.channels[c] = child1;
          sp.channels.emplace_back(child2);
        }

  return exps;
}

Solution multipleRepresentation(Solution ret) {
  memset(parent, -1, sizeof parent);
  memset(child, -1, sizeof child);

  for (int t = 0; t < ret.slots.size(); t++) {
    for (int s = 0; s < ret(t).spectrums.size(); s++) {
      Spectrum &spectrum = ret(t, s);
      int curr = 0;
      while (curr < spectrum.channels.size()) {
        if (curr + 1 < spectrum.channels.size() &&
            spectrum(curr).bandwidth == spectrum(curr + 1).bandwidth) {
          Channel merged = spectrum(curr);
          merged.bandwidth *= 2;

          for (const auto &conn : spectrum(curr + 1).connections)
            merged.connections.emplace_back(conn);

          int p = spectrum.channels.size();
          child[t][s][p][0] = curr;
          child[t][s][p][1] = curr + 1;
          parent[t][s][curr] = p;
          parent[t][s][curr + 1] = p;

          spectrum.channels.emplace_back(merged);
          curr += 2;
        } else {
          curr += 1;
        }
      }
    }
  }

  computeThroughput(ret);
  return ret;
}

Solution perturbation(Solution &sol, int kkmul) {
  int rnd = random_number(2);
  if (rnd)
    K_AddDrop(sol, kkmul);
  else
    K_RemoveAndInserts(sol, kkmul);

  computeThroughput(sol);
  fix_channels(sol);

  return sol;
}

double solve(int k, int i, int j) {
  inSolution[k][i][j] = true;
  double ret = chanOF[k][i][j];
  if (child[k][i][j][0] != -1 && child[k][i][j][1] != -1) {
    double a1 = solve(k, i, child[k][i][j][0]);
    double a2 = solve(k, i, child[k][i][j][1]);

    if (a1 + a2 > ret) {
      ret = a1 + a2;
      inSolution[k][i][j] = false;
    }
  }

  return ret;
}

void setDP(Solution &sol) {
  memset(chanOF, 0, sizeof chanOF);
  memset(inSolution, false, sizeof inSolution);

  for (int t = 0; t < sol.slots.size(); t++)
    for (int s = 0; s < sol(t).spectrums.size(); s++)
      for (int c = 0; c < sol(t, s).channels.size(); c++)
        chanOF[t][s][c] = sol(t, s, c).throughput;
}

void setDP(Solution &sol, int t, int s) {
  for (int c = 0; c < sol(t, s).channels.size(); c++) {
    inSolution[t][s][c] = false;

    chanOF[t][s][c] = sol(t, s, c).throughput;
  }
}

double calcDP(Solution &sol, int t, int s) {
  double OF = 0.0;

  for (int c = 0; c < sol(t, s).channels.size(); c++) {
    if (parent[t][s][c] == -1) {
      double ret = solve(t, s, c);
      OF += ret;
    }
  }

  return OF;
}

double calcDP(Solution &sol) {
  double OF = 0.0;

  for (int t = 0; t < sol.slots.size(); t++) {
    for (int s = 0; s < sol(t).spectrums.size(); s++) {
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        if (parent[t][s][c] == -1) {
          double ret = solve(t, s, c);
          OF += ret;
        }
      }
    }
  }

  return OF;
}

void removeAllOccurrences(Solution &sol, int id) {
  for (int t = 0; t < sol.slots.size(); t++)
    for (int s = 0; s < sol(t).spectrums.size(); s++)
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        Channel &chan = sol(t, s, c);
        for (const Connection &conn : chan.connections) {
          if (conn.id == id) {
            chan = deleteFromChannel(chan, id);
            break;
          }
        }
      }
}

void addEverywhere(Solution &sol, int id) {
  for (int t = 0; t < sol.slots.size(); t++)
    for (int s = 0; s < sol(t).spectrums.size(); s++)
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        Channel &ch = sol(t, s, c);
        ch = insertInChannel(ch, id);
      }
}

Solution local_search(Solution &mult, Solution &sol20) {
  bool improved = false;
  double highest = sol20.throughput;

  do {
    improved = false;

    vector<int> target;
    for (int i = 0; i < N; i++) target.emplace_back(i);

    for (int i : target) {
      Solution mult_clean(mult);
      removeAllOccurrences(mult_clean, i);

      Solution mult_cont(mult_clean);
      addEverywhere(mult_cont, i);

      double bestOF = highest;
      ti3 best_ch = {-1, -1, -1};

      setDP(mult_clean);
      double cleanOF = calcDP(mult);

      for (int t = 0; t < sol20.slots.size(); t++)
        for (int s = 0; s < sol20(t).spectrums.size(); s++)
          for (int c = 0; c < sol20(t, s).channels.size(); c++) {
            setDP(mult_clean);
            int c_ch = c;

            while (c_ch != -1) {
              chanOF[t][s][c_ch] = mult_cont(t, s, c_ch).throughput;
              c_ch = parent[t][s][c_ch];
            }

            double ans = calcDP(mult);
            if (compareObjectives(ans, bestOF) < 0) {
              bestOF = ans;
              best_ch = {t, s, c};
            }
            // setDP(mult_clean, t, s);
            // double AA = calcDP(mult, t, s);
            //
            // setDP(mult_clean, t, s);
            // int c_ch = c;
            // while (c_ch != -1) {
            //   chanOF[t][s][c_ch] = mult_cont(t, s, c_ch).throughput;
            //   c_ch = parent[t][s][c_ch];
            // }
            //
            // double BB = calcDP(mult, t, s);
            // double result = cleanOF - AA + BB;
            //
            // if (compareObjectives(result, bestOF) < 0) {
            //   bestOF = result;
            //   best_ch = {t, s, c};
            // }
          }

      int better = compareObjectives(bestOF, highest);
      if (better < 0 && best_ch != make_tuple(-1, -1, -1)) {
        improved = true;
        highest = bestOF;

        mult = mult_clean;
        auto &[nt, nsp, nch] = best_ch;
        while (nch != -1) {
          mult(nt, nsp).channels[nch] = insertInChannel(mult(nt, nsp, nch), i);
          nch = parent[nt][nsp][nch];
        }
      }
    }

  } while (improved);

  setDP(mult);
  calcDP(mult);

  Solution ret = reconstruct_sol(mult);
  assert(approximatelyEqual(ret.throughput, highest));
  return ret;
}

double computeConnectionThroughput(Connection &conn, int bandwidth) {
  if (bandwidth == 0) return 0.0;

  int mcs = -1;
  int maxDataRate = 11;

  if (approximatelyEqual(conn.interference, 0.0, EPS)) {
    mcs = maxDataRate;
    conn.throughput = DR[mcs][bwIdx(bandwidth)];
  } else {
    conn.SINR = AFF[conn.id][conn.id] / (conn.interference + noise);

    mcs = 11;
    while (mcs >= 0 && conn.SINR < SINR[mcs][bwIdx(bandwidth)]) mcs--;

    if (mcs < 0)
      conn.throughput = 0.0;
    else
      conn.throughput = DR[mcs][bwIdx(bandwidth)];
  }

  return conn.throughput;
}

Channel insertInChannel(Channel newChannel,
                        int idConn) {  // TODO: remover macro
  Connection conn(idConn, 0.0, 0.0, DM[idConn][idConn]);

  for (auto &connection : newChannel.connections) {
    connection.interference += AFF[connection.id][idConn];
    conn.interference += AFF[idConn][connection.id];
  }

  newChannel.connections.emplace_back(conn);
  newChannel.throughput = 0.0;
  newChannel.violation = 0.0;
  for (auto &connection : newChannel.connections) {
    computeConnectionThroughput(connection, newChannel.bandwidth);
    newChannel.throughput += connection.throughput;
  }

  return newChannel;
}

Channel deleteFromChannel(const Channel &channel, int idConn) {
  Channel newChannel(channel.bandwidth);

  for (const Connection &conn : channel.connections)
    if (conn.id != idConn) newChannel.connections.emplace_back(conn);

  newChannel.throughput = 0.0;
  newChannel.violation = 0.0;
  for (Connection &conn : newChannel.connections) {
    conn.interference -= AFF[conn.id][idConn];

    computeConnectionThroughput(conn, newChannel.bandwidth);
    newChannel.throughput += conn.throughput;
#ifdef MDVRBSP
    // TODO: conectar connection.violation aqui
    if (of_opt == SUMVIO)
      newChannel.violation += max(0.0, gma[conn.id] - conn.throughput);
    else if (of_opt == MINMAX)
      newChannel.violation =
          max(newChannel.violation, gma[conn.id] - conn.throughput);
#endif
  }

  return newChannel;
}

bool reinsert(Solution &sol, Connection conn, ti3 from, ti3 to, bool force) {
  if (from == to) return false;

  const auto &[ft, fs, fc] = from;
  const auto &[tt, ts, tc] = to;

  Channel &new_chan = sol(tt, ts, tc);
  Channel &old_chan = sol(ft, fs, fc);

  Channel old_new_chan = deleteFromChannel(old_chan, conn.id);
  Channel copy_new_chan = insertInChannel(new_chan, conn.id);

  double newObjective = sol.throughput - old_chan.throughput +
                        old_new_chan.throughput - new_chan.throughput +
                        copy_new_chan.throughput;

  bool improved = false;
  if (newObjective > sol.throughput) improved = true;

  if ((newObjective > sol.throughput) || force) {
    old_chan = old_new_chan;
    new_chan = copy_new_chan;
    sol.throughput = newObjective;
  }

  return improved;
}

double computeThroughput(Solution &curr, bool force) {
  double OF = 0.0;

  for (int t = 0; t < curr.slots.size(); t++) {
    for (int s = 0; s < curr(t).spectrums.size(); s++) {
      for (int c = 0; c < curr(t, s).channels.size(); c++) {
        if (make_tuple(t, s, c) == zeroChannel) continue;

        auto &chan = curr(t, s, c);
        chan.throughput = 0.0;
        for (Connection &u : chan.connections) {
          u.interference = 0.0;
          u.throughput = 0.0;
          for (Connection &v : chan.connections)
            u.interference += u.id == v.id ? 0.0 : AFF[u.id][v.id];

          chan.throughput += computeConnectionThroughput(u, chan.bandwidth);
        }
        OF += chan.throughput;
      }
    }
  }

  curr.throughput = OF;
  return OF;
}

int bwIdx(int bw) {
  if (bw == 40) {
    return 1;
  } else if (bw == 80) {
    return 2;
  } else if (bw == 160) {
    return 3;
  }
  return 0;
}

Solution reconstruct_sol(Solution &curr) {
  Solution ret;
  ret.slots.resize(curr.slots.size());

  for (int t = 0; t < curr.slots.size(); t++)
    ret(t).spectrums.resize(curr(t).spectrums.size(),
                            Spectrum(0.0, 0.0, vector<Channel>()));

  for (int t = 0; t < curr.slots.size(); t++) {
    for (int s = 0; s < curr(t).spectrums.size(); s++) {
      for (int c = 0; c < curr(t, s).channels.size(); c++) {
        if (parent[t][s][c] == -1) {
          recoverSolution(t, s, c, false);
        }
      }
    }
  }

  for (int t = 0; t < MAX_SLOTS; t++) {
    for (int s = 0; s < MAX_SPECTRUM; s++) {
      for (int c = 0; c < MAX_CHANNELS; c++) {
        if (inSolution[t][s][c]) {
          ret(t, s).channels.emplace_back(curr(t, s, c));
        }
      }
    }
  }

  computeThroughput(ret);

  return ret;
}

void recoverSolution(int k, int i, int j, bool clean) {
  if (clean) inSolution[k][i][j] = false;

  if (inSolution[k][i][j]) {
    clean = true;
  }

  if (child[k][i][j][0] != -1 && child[k][i][j][1] != -1) {
    recoverSolution(k, i, child[k][i][j][0], clean);
    recoverSolution(k, i, child[k][i][j][1], clean);
  }
}

void convertTableToMW(const vector<vector<double>> &_SINR,
                      vector<vector<double>> &__SINR) {
  double result, b;
  for (int i = 0; i < __SINR.size(); i++) {
    for (int j = 0; j < __SINR[i].size(); j++) {
      if (_SINR[i][j] != 0) {
        b = _SINR[i][j] / 10.0;  // dBm divided by 10
        result = pow(10.0, b);   // Convert DBM to mW

        __SINR[i][j] = result;
      } else {
        __SINR[i][j] = 0;
      }
    }
  }
}

void distanceAndInterference() {
  for (int i = 0; i < N; i++) {
    double sender_i_x = senders[i][0];
    double sender_i_y = senders[i][1];

    for (int j = 0; j < N; j++) {
      double receiver_j_x = receivers[j][0];
      double receiver_j_y = receivers[j][1];

      double dist =
          distance(sender_i_x, sender_i_y, receiver_j_x, receiver_j_y);
      DM[i][j] = dist;
    }
  }
}

double convertDBMToMW(double _value) {
  double result = 0.0, b;

  b = _value / 10.0;      // dBm dividido por 10
  result = pow(10.0, b);  // Converte de DBm para mW

  return result;
}

double gammaToBeta(double gmm, int bw) {
  int last = -1;

  for (int i = 0; i < 12; i++) {
    if (DR[i][bw] >= gmm) {
      last = i;
      break;
    }
  }

  if (last != -1) {
    return SINR[last][bw];
  }

  return SINR[11][3] + 1.0;
}

void computeChannelsThroughput(vector<Channel> &channels) {
  for (Channel &channel : channels) {
    double of = 0.0;
    for (Connection &conn : channel.connections) {
      of += computeConnectionThroughput(conn, channel.bandwidth);
    }
    channel.throughput = of;
  }
}

optional<vector<Channel>> split(Channel toSplit, vector<double> *GMA) {
  int newBw = toSplit.bandwidth / 2;
  vector<Channel> ret;
  ret.emplace_back(Channel(0.0, 0.0, 0.0, newBw, vector<Connection>()));
  ret.emplace_back(Channel(0.0, 0.0, 0.0, newBw, vector<Connection>()));

  vector<Connection> connToInsert = toSplit.connections;
  sort(connToInsert.rbegin(), connToInsert.rend());

  double bestThroughputsoFar = 0.0;
  for (int i = 0; i < connToInsert.size(); i++) {
    double bestThroughputIteration = 0.0, testando = 0.0;
    Channel bestChannel(0.0, 0.0, 0.0, newBw, vector<Connection>());
    int bestIdx = -1;

    for (int c = 0; c < ret.size(); c++) {
      Channel inserted = insertInChannel(ret[c], connToInsert[i].id);
      double resultThroughput = inserted.throughput;

      double g_throughput =
          (bestThroughputsoFar - ret[c].throughput) + resultThroughput;
      testando = max(testando, g_throughput);
      bool one = g_throughput > bestThroughputIteration;
      bool two = approximatelyEqual(g_throughput, bestThroughputIteration) &&
                 (inserted.connections.size() < bestChannel.connections.size());

      if (one || two) {
        bestThroughputIteration = g_throughput;
        bestChannel = inserted;
        bestIdx = c;
      }
    }

    assert(bestIdx >= 0);
    ret[bestIdx] = bestChannel;
    bestThroughputsoFar = bestThroughputIteration;
  }
  computeChannelsThroughput(ret);

  assert(approximatelyEqual(ret[0].throughput + ret[1].throughput,
                            bestThroughputsoFar));

  if (GMA != nullptr)
    for (const auto &ch : ret)
      for (const auto &conn : ch.connections)
        if (definitelyLessThan(conn.throughput, (*GMA)[conn.id]))
          return nullopt;

  return ret;
}

optional<vector<Channel>> split2(Channel toSplit,
                                 vector<double> *GMA = nullptr) {
  if (toSplit.bandwidth <= 20) return nullopt;

  int new_bw = toSplit.bandwidth / 2;
  vector<Connection> conns = toSplit.connections;
  vector<Channel> ret(2, Channel(0.0, 0.0, 0.0, new_bw, {}));

  for (const auto &con : conns) {
    Channel best(0.0, 0.0, 0.0, new_bw, {});
    int it_best = -1;

    for (int c = 0; c < ret.size(); ++c) {
      Channel ans = insertInChannel(ret[c], con.id);

      if (definitelyGreaterThan(ans.throughput, best.throughput)) {
        best = ans;
        it_best = c;
      }
    }

    assert(it_best >= 0);
    ret[it_best] = best;
  }

  if (GMA != nullptr)
    for (const auto &ch : ret)
      for (const auto &conn : ch.connections)
        if (definitelyLessThan(conn.throughput, (*GMA)[conn.id]))
          return nullopt;

  return ret;
}

Solution CH_VRBSP() {
  auto compare_sol = [](Solution *a, Solution *b) {
    return definitelyLessThan(a->throughput, b->throughput);
  };

  Solution ret(init_conf, 0.0, true);
  vector<int> links;
  for (int i = 0; i < N; i++) links.emplace_back(i);

  shuffle(links.begin(), links.end(), gen);

  for (int i = 0; i < N; ++i) {
    Solution tmp = ret;

    for (int s = 0; s < ret(0).spectrums.size(); ++s)
      for (int c = 0; c < ret(0, s).channels.size(); ++c) {
        Channel c1 = insertInChannel(ret(0, s, c), links[i]);

        Solution s1 = ret;
        s1(0, s, c) = c1;

        optional<vector<Channel>> opt_channels = *split2(ret(0, s, c));

        if (!opt_channels.value().empty()) {
          vector<Channel> &channels = *opt_channels;
          Channel c2 = insertInChannel(channels[0], links[i]);
          Channel c3 = insertInChannel(channels[1], links[i]);

          Solution s2 = ret;
          s2(0, s, c) = c2;
          s2(0, s).channels.emplace_back(channels[1]);

          Solution s3 = ret;
          s3(0, s, c) = c3;
          s3(0, s).channels.emplace_back(channels[0]);

          Solution *aux[] = {&tmp, &s1, &s2, &s3};
          for (int s = 1; s < 4; ++s) computeThroughput(*aux[s]);

          tmp = **max_element(begin(aux), end(aux), compare_sol);
          continue;
        }

        computeThroughput(s1);
        if (definitelyGreaterThan(s1.throughput, tmp.throughput)) tmp = s1;
      }

    ret = tmp;

    set<int> seen;
    for (const auto &ts : ret.slots)
      for (const auto &sp : ts.spectrums)
        for (const auto &ch : sp.channels)
          for (const auto &conn : ch.connections) {
            assert(seen.find(conn.id) == seen.end());
            seen.insert(conn.id);
          }
  }

  return ret;
}

void rmvUnusedTS(Solution &sol) {
  set<int> rmv;
  for (int t = 0; t < int(sol.slots.size()); t++) {
    bool hasConn = false;
    for (int s = 0; s < int(sol(t).spectrums.size()); s++)
      for (int c = 0; c < int(sol(t, s).channels.size()); c++)
        if (make_tuple(t, s, c) == zeroChannel)
          continue;
        else
          hasConn |= !sol(t, s, c).connections.empty();

    if (!hasConn) rmv.insert(t);
  }

  Solution aux;
  for (int t = 0; t < int(sol.slots.size()); t++)
    if (rmv.count(t) == 0) aux.slots.emplace_back(sol(t));

  aux.violation = sol.violation;
  aux.throughput = sol.throughput;

  if (!rmv.empty()) printf("rmvUnusedTS removed %lu ts\n", rmv.size());
}

void print_solution_to_gurobi(const Solution &sol, FILE *file) {
  // assert(sol.slots[0].spectrums[3].channels[0].connections.empty());
  // if (file == nullptr) file = fopen("./solution_gurobi.mst", "w");
  //
  // vector<set<int>> aux(4, set<int>());
  // for (int j = 0; j < 45; j++) aux[bwIdx2(j)].insert(j);
  //
  // int cnt = 0;
  // vector<pair<vector<int>, int>> S;
  // for (auto &ts : sol.slots)
  //   for (auto &sp : ts.spectrums)
  //     for (auto &ch : sp.channels) {
  //       if (ch.bandwidth == 0) continue;
  //
  //       vector<int> conn_v;
  //       for (auto &conn : ch.connections) conn_v.emplace_back(conn.id);
  //
  //       cnt += int(conn_v.size());
  //       S.emplace_back(make_pair(conn_v, ch.bandwidth));
  //     }
  //
  // sort(S.rbegin(), S.rend(),
  //      [](const pair<vector<int>, int> &a, const pair<vector<int>, int> &b) {
  //        return a.second < b.second;
  //      });
  //
  // int bw_cnt[] = {0, 0, 0, 0};
  //
  // for (auto s : S) ++bw_cnt[bwIdx(s.second)];
  //
  // assert(cnt == n_connections);
  //
  // // fprintf(file, "# %lu\n", sol.slots.size());
  //
  // vector<vector<set<int>>> T;
  // T.assign(sol.slots.size(), aux);
  //
  // for (const auto &s : S) {
  //   int bw = s.second;
  //   int bw_idx = bwIdx(bw);
  //
  //   bool some = false;
  //   for (int t_ = 0; t_ < int(T.size()); ++t_) {
  //     bool found = false;
  //     if (!T[t_][bw_idx].empty()) {
  //       int ch_idx = *T[t_][bw_idx].begin();
  //
  //       for (const auto &c : s.first)
  //         fprintf(file, "%d %d %d 1\n", c, ch_idx, t_);
  //
  //       for (int j = 0; j < 45; ++j)
  //         if (overlap[ch_idx][j] && T[t_][bwIdx2(j)].count(j))
  //           T[t_][bwIdx2(j)].erase(j);
  //
  //       --bw_cnt[bw_idx];
  //       found = true;
  //       some = true;
  //     }
  //
  //     if (found) break;
  //   }
  //
  //   assert(some);
  // }
  //
  // fclose(file);
  // file = nullptr;
}

void initTimeSlot() {
  for (Spectrum &sp : init_conf) {
    while (sp.maxFrequency - sp.usedFrequency > 0) {
      int bw = 160;
      while (bw > (sp.maxFrequency - sp.usedFrequency) && bw > 20) {
        bw /= 2;
      }

      if (bw <= (sp.maxFrequency - sp.usedFrequency)) {
        sp.usedFrequency += bw;
        sp.channels.emplace_back(0.0, 0.0, 0.0, bw, vector<Connection>());
      }
    }
    assert(sp.maxFrequency - sp.usedFrequency >= 0);
  }

  // vector<Channel> aux;
  // aux.emplace_back(Channel());
  // init_conf.emplace_back(0, 0, aux);
}

#ifndef MDVRBSP
void read_data() {
#else
void read_data(vector<double> &GMA) {
#endif
  scanf("%d %lf %lf %lf %d", &N, &alfa, &noise, &powerSender, &n_spectrums);

  for (int i = 0; i < n_spectrums; i++) {
    int s;
    scanf("%d", &s);
    spectrum_size.emplace_back(s);
    init_conf.emplace_back(s, 0, vector<Channel>());
  }

  if (noise != 0) {
    noise = convertDBMToMW(noise);
  }

  for (int i = 0; i < N; i++) {
    double a, b;
    scanf("%lf %lf", &a, &b);
    receivers[i][0] = a;
    receivers[i][1] = b;
  }

  for (int i = 0; i < N; i++) {
    double a, b;
    scanf("%lf %lf", &a, &b);
    senders[i][0] = a;
    senders[i][1] = b;
  }

  for (int i = 0; i < N; i++) {
    double bt;
    scanf("%lf", &bt);
#ifdef MDVRBSP
    GMA.emplace_back(bt);
#endif
  }

  DR.assign(12, vector<double>(4, 0));
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 4; j++) {
      double a;
      scanf("%lf", &a);
      DR[i][j] = a;
    }
  }

  SINR.assign(12, vector<double>(4, 0));
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 4; j++) {
      double a;
      scanf("%lf", &a);
      SINR[i][j] = a;
    }
  }

  convertTableToMW(SINR, SINR);
  distanceAndInterference();
  initTimeSlot();

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double value = powerSender / pow(DM[i][j], alfa);
      AFF[i][j] = value;
    }
  }
}

void print_objective_to_file(const Solution &sol, FILE *file) {
  if (file == nullptr) {
    file = fopen("./objectives.txt", "a");
  }

#ifdef MDVRBSP
  fprintf(file, "%lf\n", sol.violation);
#else
  fprintf(file, "%lf\n", sol.throughput);
#endif

  fclose(file);
}

int compareObjectives(const int lhs, const int rhs) {
  if (definitelyGreaterThan(lhs, rhs)) return -1;
  return essentiallyEqual(lhs, rhs) ? 0 : 1;
}

int compareObjectives(const Solution &lhs, const Solution &rhs) {
  return compareObjectives(lhs.throughput, rhs.throughput);
}

Connection::Connection(int id, double throughput, double interference,
                       double distanceSR)
    : id(id),
      throughput(throughput),
      interference(interference),
      distanceSR(distanceSR) {
  violation = 0.0;
  SINR = 0.0;
}

Connection::Connection(int id) : id(id) {
  throughput = 0.0;
  interference = 0.0;
  SINR = 0.0;
  violation = 0.0;
  distanceSR = DM[id][id];
}

bool Connection::operator<(const Connection &other) const {
  return distanceSR < other.distanceSR;
}

bool Connection::operator>(const Connection &other) const {
  return !operator<(other);
}

bool Channel::operator<(const Channel &other) const {
  return bandwidth < other.bandwidth;
}

bool Channel::operator>(const Channel &other) const {
  return !operator<(other);
}

Channel::Channel(double throughput, double interference, double violation,
                 int bandwidth, const vector<Connection> connections)
    : throughput(throughput),
      interference(interference),
      violation(violation),
      bandwidth(bandwidth),
      connections(connections) {}

Channel::Channel(int bandwidth)
    : Channel(0.0, 0.0, 0.0, bandwidth, vector<Connection>()) {}

Channel::Channel() : Channel(0.0) {}

Connection &Channel::operator()(int l) { return connections[l]; }

Spectrum::Spectrum(int maxFrequency, int usedFrequency,
                   const vector<Channel> channels)
    : maxFrequency(maxFrequency),
      usedFrequency(usedFrequency),
      channels(channels) {
  interference = 0.0;
}

Spectrum::Spectrum() {
  maxFrequency = 0;
  usedFrequency = 0;
  interference = 0.0;
  channels = vector<Channel>();
}

Channel &Spectrum::operator()(int c) { return channels[c]; }

Connection &Spectrum::operator()(int c, int l) {
  return channels[c].connections[l];
}

TimeSlot::TimeSlot(const vector<Spectrum> sp) : TimeSlot() { spectrums = sp; }

TimeSlot::TimeSlot() {
  interference = 0.0;
  throughput = 0.0;
  spectrums = vector<Spectrum>();
}

Spectrum &TimeSlot::operator()(int s) { return spectrums[s]; }

Channel &TimeSlot::operator()(int s, int c) { return spectrums[s].channels[c]; }

Connection &TimeSlot::operator()(int s, int c, int l) {
  return spectrums[s].channels[c].connections[l];
}

Solution::Solution(const vector<Spectrum> sp, double tot, bool flag)
    : throughput(tot) {
  slots.emplace_back(sp);
  violation = 0.0;
}

Solution::Solution(const vector<TimeSlot> &ts) : slots(ts) {
  throughput = 0.0;
  violation = 0.0;
}

Solution::Solution() {
  slots = vector<TimeSlot>();
  throughput = 0.0;
  violation = 0.0;
}

#ifdef MDVRBSP
bool Solution::operator<(const Solution &o1) const {
  return violation < o1.violation;
}

bool Solution::operator>(const Solution &o1) const {
  return !Solution::operator<(o1);
}
#else
bool Solution::operator<(const Solution &o1) const {
  return throughput < o1.throughput;
}

bool Solution::operator>(const Solution &o1) const {
  return !Solution::operator<(o1);
}
#endif

TimeSlot &Solution::operator()(int t) { return slots[t]; }

Spectrum &Solution::operator()(int t, int s) { return slots[t].spectrums[s]; }

Channel &Solution::operator()(int t, int s, int c) {
  return slots[t].spectrums[s].channels[c];
}

Connection &Solution::operator()(int t, int s, int c, int l) {
  return slots[t].spectrums[s].channels[c].connections[l];
}
