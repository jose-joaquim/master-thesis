#ifndef COMMON_H
#define COMMON_H

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <thread>
#include <tuple>
#include <vector>

const int MAX_CONN = 2048;
const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-9;

using namespace std;
using namespace std::chrono;

enum CH_PARAMS { RANDOM, GREEDY, MIN_GAMMA, MINIMUM_AFF };
#ifdef MDVRBSP
enum OF_TYPE { MINMAX, SUMVIO };
#endif

using ii = pair<int, int>;
using dd = pair<double, double>;
using ti3 = tuple<int, int, int>;

static int n_connections, n_spectrums;
static double alfa, noise, powerSender;
static double receivers[MAX_CONN][2], senders[MAX_CONN][2];
static double affectance[MAX_CONN][MAX_CONN];
static double distanceMatrix[MAX_CONN][MAX_CONN];
static vector<vector<double>> dataRates, SINR;
static vector<int> spectrum_size;

static const int MAX_SPECTRUM = 5;
static const int MAX_CHANNELS = 45;
static const int MAX_SLOTS = 2048;
// static double startTime;
static double maximumTime;

static high_resolution_clock::time_point start;

static int parent[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
static int child[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS][2];
static double chanOF[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
static bool inSolution[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
static CH_PARAMS ch_opt;

#ifdef MDVRBSP
static OF_TYPE of_opt;
static vector<double> gma;
#endif

random_device rd;
auto whatever = default_random_engine{rd()};
uniform_real_distribution<int> dist(0, MAX_CONN);
static ti3 zeroChannel = {0, 3, 0};

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

struct Connection {
  int id;
  double throughput;
  double interference;
  double SINR;
  double distanceSR;
  double violation;

  Connection(int id, double throughput, double interference, double distanceSR)
      : id(id),
        throughput(throughput),
        interference(interference),
        distanceSR(distanceSR) {
    violation = 0.0;
    SINR = 0.0;
  }

  Connection(int id) : id(id) {
    throughput = 0.0;
    interference = 0.0;
    SINR = 0.0;
    violation = 0.0;
    distanceSR = distanceMatrix[id][id];
  }

  bool operator<(const Connection &other) const {
    return distanceSR < other.distanceSR;
  }

  bool operator>(const Connection &other) const { return !operator<(other); }
};

struct Channel {
  double throughput;
  double interference;
  double violation;
  int bandwidth;
  vector<Connection> connections;

  bool operator<(const Channel &other) const {
    return bandwidth < other.bandwidth;
  }

  bool operator>(const Channel &other) const { return !operator<(other); }

  Channel(double throughput, double interference, double violation,
          int bandwidth, const vector<Connection> connections)
      : throughput(throughput),
        interference(interference),
        violation(violation),
        bandwidth(bandwidth),
        connections(connections) {}

  Channel(int bandwidth)
      : Channel(0.0, 0.0, 0.0, bandwidth, vector<Connection>()) {}

  Channel() : Channel(0.0) {}

  Connection &operator()(int l) { return connections[l]; }
};

struct Spectrum {
  int maxFrequency;
  int usedFrequency;
  double interference;
  vector<Channel> channels;

  Spectrum(int maxFrequency, int usedFrequency, const vector<Channel> channels)
      : maxFrequency(maxFrequency),
        usedFrequency(usedFrequency),
        channels(channels) {
    interference = 0.0;
  }

  Spectrum() {
    maxFrequency = 0;
    usedFrequency = 0;
    interference = 0.0;
    channels = vector<Channel>();
  }

  Channel &operator()(int c) { return channels[c]; }
  Connection &operator()(int c, int l) { return channels[c].connections[l]; }
};

struct TimeSlot {
  vector<Spectrum> spectrums;
  double interference;
  double throughput;

  TimeSlot(const vector<Spectrum> sp) : TimeSlot() { spectrums = sp; }
  TimeSlot() {
    interference = 0.0;
    throughput = 0.0;
    spectrums = vector<Spectrum>();
  }

  Spectrum &operator()(int s) { return spectrums[s]; }
  Channel &operator()(int s, int c) { return spectrums[s].channels[c]; }
  Connection &operator()(int s, int c, int l) {
    return spectrums[s].channels[c].connections[l];
  }
};

class Solution {
 public:
  vector<TimeSlot> slots;
  double throughput;
  double violation;

  Solution(const vector<Spectrum> sp, double tot, bool flag) : throughput(tot) {
    slots.emplace_back(sp);
    violation = 0.0;
  }

  Solution(const vector<TimeSlot> &ts) : slots(ts) {
    throughput = 0.0;
    violation = 0.0;
  }

  Solution() {
    slots = vector<TimeSlot>();
    throughput = 0.0;
    violation = 0.0;
  }

#ifdef MDVRBSP
  bool operator<(const Solution &o1) const { return violation < o1.violation; }

  bool operator>(const Solution &o1) const { return !operator<(o1); }
#else
  bool operator<(const Solution &o1) const {
    return throughput < o1.throughput;
  }

  bool operator>(const Solution &o1) const { return !operator<(o1); }
#endif

  TimeSlot &operator()(int t) { return slots[t]; }
  Spectrum &operator()(int t, int s) { return slots[t].spectrums[s]; }
  Channel &operator()(int t, int s, int c) {
    return slots[t].spectrums[s].channels[c];
  }
  Connection &operator()(int t, int s, int c, int l) {
    return slots[t].spectrums[s].channels[c].connections[l];
  }
};

void init(void);

Solution vns(Solution, string);

Solution vns(string);

Solution constructive_heuristic(void);

Solution constructive_heuristic(string);

bool is_feasible(const Solution &, bool);

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

bool approximatelyEqual(double a, double b, double epsilon = EPS) {
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentiallyEqual(double a, double b, double epsilon = EPS) {
  return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyGreaterThan(double a, double b, double epsilon = EPS) {
  return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(double a, double b, double epsilon = EPS) {
  return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

void checkRepeat(const Channel &chan) {
  set<int> seen;
  for (auto x : chan.connections) {
    if (seen.count(x.id) > 0) {
      puts("opa");
    }

    seen.insert(x.id);
  }
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

bool stop(void) {
  return duration_cast<duration<double>>(high_resolution_clock::now() - start)
             .count() >= maximumTime;
}

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
  return hypot((X_si - X_ri), (Y_si - Y_ri));
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
  for (int i = 0; i < n_connections; i++) {
    double sender_i_x = senders[i][0];
    double sender_i_y = senders[i][1];

    for (int j = 0; j < n_connections; j++) {
      double receiver_j_x = receivers[j][0];
      double receiver_j_y = receivers[j][1];

      double dist =
          distance(sender_i_x, sender_i_y, receiver_j_x, receiver_j_y);
      distanceMatrix[i][j] = dist;
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
    if (dataRates[i][bw] >= gmm) {
      last = i;
      break;
    }
  }

  if (last != -1) {
    return SINR[last][bw];
  }

  return SINR[11][3] + 1.0;
}

int numberToBw(int c) {
  if (c <= 24)
    return 20;
  else if (c <= 36)
    return 40;
  else if (c <= 42)
    return 80;
  else
    return 160;
}

int cToBIdx(int c) {
  if (c <= 24)
    return 0;
  else if (c <= 36)
    return 1;
  else if (c <= 42)
    return 2;
  else
    return 3;
}

int bToIdx(int b) {
  if (b == 20)
    return 0;
  else if (b == 40)
    return 1;
  else if (b == 80)
    return 2;

  return 3;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static vector<Spectrum> init_conf;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  scanf("%d %lf %lf %lf %d", &n_connections, &alfa, &noise, &powerSender,
        &n_spectrums);

  for (int i = 0; i < n_spectrums; i++) {
    int s;
    scanf("%d", &s);
    spectrum_size.emplace_back(s);
    init_conf.emplace_back(s, 0, vector<Channel>());
  }

  if (noise != 0) {
    noise = convertDBMToMW(noise);
  }

  for (int i = 0; i < n_connections; i++) {
    double a, b;
    scanf("%lf %lf", &a, &b);
    receivers[i][0] = a;
    receivers[i][1] = b;
  }

  for (int i = 0; i < n_connections; i++) {
    double a, b;
    scanf("%lf %lf", &a, &b);
    senders[i][0] = a;
    senders[i][1] = b;
  }

  for (int i = 0; i < n_connections; i++) {
    double bt;
    scanf("%lf", &bt);
#ifdef MDVRBSP
    GMA.emplace_back(bt);
#endif
  }

  dataRates.assign(12, vector<double>(4, 0));
  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 4; j++) {
      double a;
      scanf("%lf", &a);
      dataRates[i][j] = a;
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

  for (int i = 0; i < n_connections; i++) {
    for (int j = 0; j < n_connections; j++) {
      double value = powerSender / pow(distanceMatrix[i][j], alfa);
      affectance[i][j] = value;
    }
  }
}

void print_objective_to_file(const Solution &sol, FILE *file = nullptr) {
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

int bwIdx2(int x) {
  if (x < 25) return 0;
  if (x < 37) return 1;
  if (x < 43) return 2;

  return 3;
}

int count_conn(const Solution &, bool op = false);

void print_solution_to_gurobi(const Solution &sol, FILE *file = nullptr) {
  assert(sol.slots[0].spectrums[3].channels[0].connections.empty());
  if (file == nullptr) file = fopen("./solution_gurobi.mst", "w");

  vector<set<int>> aux(4, set<int>());
  for (int j = 0; j < 45; j++) aux[bwIdx2(j)].insert(j);

  int cnt = 0;
  vector<pair<vector<int>, int>> S;
  for (auto &ts : sol.slots)
    for (auto &sp : ts.spectrums)
      for (auto &ch : sp.channels) {
        if (ch.bandwidth == 0) continue;

        vector<int> conn_v;
        for (auto &conn : ch.connections) conn_v.emplace_back(conn.id);

        cnt += int(conn_v.size());
        S.emplace_back(make_pair(conn_v, ch.bandwidth));
      }

  sort(S.rbegin(), S.rend(),
       [](const pair<vector<int>, int> &a, const pair<vector<int>, int> &b) {
         return a.second < b.second;
       });

  int bw_cnt[] = {0, 0, 0, 0};

  for (auto s : S) ++bw_cnt[bwIdx(s.second)];

  assert(cnt == n_connections);

  // fprintf(file, "# %lu\n", sol.slots.size());

  vector<vector<set<int>>> T;
  T.assign(sol.slots.size(), aux);

  // for (int t = 0; t < int(sol.slots.size()); t++) fprintf(file, "t[%d] 1\n",
  // t);

  for (const auto &s : S) {
    int bw = s.second;
    int bw_idx = bwIdx(bw);

    bool some = false;
    for (int t_ = 0; t_ < int(T.size()); ++t_) {
      bool found = false;
      if (!T[t_][bw_idx].empty()) {
        int ch_idx = *T[t_][bw_idx].begin();

        for (const auto &c : s.first)
          fprintf(file, "%d %d %d 1\n", c, ch_idx, t_);

        for (int j = 0; j < 45; ++j)
          if (overlap[ch_idx][j] && T[t_][bwIdx2(j)].count(j))
            T[t_][bwIdx2(j)].erase(j);

        --bw_cnt[bw_idx];
        found = true;
        some = true;
      }

      if (found) break;
    }

    assert(some);
  }

  fclose(file);
  file = nullptr;
  assert(count_conn(sol));
}

void print_solution_to_file(Solution &sol, FILE *file = nullptr) {
  if (file == nullptr) file = fopen("./solution.txt", "w");

  fprintf(file, "\n%lf %d\n", sol.throughput, int(sol.slots.size()));
  for (int t = 0; t < sol.slots.size(); t++) {
    fprintf(file, "%d\n", int(sol(t).spectrums.size()));
    for (int s = 0; s < sol(t).spectrums.size(); s++) {
      fprintf(file, "%d\n", int(sol(t, s).channels.size()));
      for (int c = 0; c < int(sol(t, s).channels.size()); c++) {
        fprintf(file, "%d %d", sol(t, s, c).bandwidth,
                int(sol(t, s, c).connections.size()));
        for (int i = 0; i < sol(t, s, c).connections.size(); i++) {
          fprintf(file, " %d", sol(t, s, c).connections[i].id);
        }
        fprintf(file, "\n");
      }
    }
  }

  fclose(file);
}

double computeConnectionThroughput(Connection &conn, int bandwidth) {
  if (bandwidth == 0) return 0.0;

  int mcs = -1;
  int maxDataRate = 11;

  if (approximatelyEqual(conn.interference, 0.0, EPS)) {
    mcs = maxDataRate;
    conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
  } else {
    // double conn_SINR = (powerSender / pow(conn.distanceSR, alfa)) /
    // (conn.interference + noise); conn.SINR = conn_SINR;

    conn.SINR = affectance[conn.id][conn.id] / (conn.interference + noise);

    mcs = 11;
    while (mcs >= 0 && conn.SINR < SINR[mcs][bwIdx(bandwidth)]) mcs--;

    if (mcs < 0)
      conn.throughput = 0.0;
    else
      conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
  }

  return conn.throughput;
}

Channel insertInChannel(Channel newChannel,
                        int idConn) {  // TODO: remover macro
  Connection conn(idConn, 0.0, 0.0, distanceMatrix[idConn][idConn]);

  for (auto &connection : newChannel.connections) {
    connection.interference += affectance[connection.id][idConn];
    conn.interference += affectance[idConn][connection.id];
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
    conn.interference -= affectance[conn.id][idConn];

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

bool reinsert(Solution &sol, Connection conn, ti3 from, ti3 to,
              bool force = false) {
  if (from == to) return false;

  const auto &[ft, fs, fc] = from;
  const auto &[tt, ts, tc] = to;

  Channel &new_chan = sol(tt, ts, tc);
  Channel &old_chan = sol(ft, fs, fc);

  checkRepeat(new_chan);
  checkRepeat(old_chan);

  Channel old_new_chan = deleteFromChannel(old_chan, conn.id);
  Channel copy_new_chan = insertInChannel(new_chan, conn.id);

  checkRepeat(old_new_chan);
  checkRepeat(copy_new_chan);

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

  checkRepeat(new_chan);

  return improved;
}

double computeThroughput(Solution &curr, bool force = false) {
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
            u.interference += u.id == v.id ? 0.0 : affectance[u.id][v.id];

          chan.throughput += computeConnectionThroughput(u, chan.bandwidth);
        }
        OF += chan.throughput;
      }
    }
  }

  curr.throughput = OF;
  return OF;
}

void print_solution(Solution &sol, bool ok = true) {
  printf("sol w/ %lu time-slots, throughput %.2lf, violation %.2lf\n",
         sol.slots.size(), sol.throughput, sol.violation);
  for (int t = 0; t < sol.slots.size(); t++) {
    printf("TS: %d:\n{\n", t);
    for (int s = 0; s < sol(t).spectrums.size() - 1; s++) {
      printf("Spec: %d, ", s);
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        printf("CH: %d (bw: %d)\nconns:\n", c, sol(t, s, c).bandwidth);
        for (const Connection &conn : sol(t, s, c).connections) {
          printf("{%d, %lf, %lf, %lf}\n", conn.id, conn.throughput,
                 conn.interference, conn.SINR);
        }
        puts("");
      }
    }
    printf("}\n");
  }

  puts("unscheduled connections:");
  for (const Connection &c :
       sol.slots[0].spectrums[3].channels[0].connections) {
    printf(" %d", c.id);
  }
  puts("");
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

int count_conn(const Solution &sol, bool op) {
  int ret = 0;
#ifdef MDVRBSP
  if (!sol.slots[0].spectrums[3].channels[0].connections.empty()) {
    for (auto &c : sol.slots[0].spectrums[3].channels[0].connections)
      printf("%d ", c.id);
    puts("");
  }

  assert(sol.slots[0].spectrums[3].channels[0].connections.empty());
#endif

  for (auto &x : sol.slots)
    for (auto &y : x.spectrums)
      for (auto &z : y.channels)
        // if (z.bandwidth != 0)
        ret += z.connections.size();

  bool fail = op ? ret < n_connections : ret != n_connections;
  if (fail) {
    printf(" vixe %d %d\n", ret, n_connections);

    set<int> ids;
    for (auto &x : sol.slots)
      for (auto &y : x.spectrums)
        for (auto &z : y.channels)
          for (auto &c : z.connections) {
            ids.insert(c.id);
            printf(" %d", c.id);
          }

    puts("");
    for (int kk = 0; kk < n_connections; kk++) {
      if (ids.count(kk) == 0) printf("faltou %d\n", kk);
    }

    puts("");
    assert(false);
  }
  return ret;
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
  //     computeThroughput(sol);
  //
  // #ifdef MDVRBSP
  //     computeViolation(sol);
  // #endif
  return improved;
}

Solution convertTo20MhzSol(Solution exps) {
  for (TimeSlot &ts : exps.slots)
    for (Spectrum &sp : ts.spectrums)
      sort(sp.channels.rbegin(), sp.channels.rend());

  while (!allChannels20MHz(exps)) {
    for (TimeSlot &ts : exps.slots) {
      for (Spectrum &sp : ts.spectrums) {
        for (int c = 0; c < sp.channels.size(); c++) {
          if (sp.channels[c].bandwidth <= 20) continue;

          int newBw = sp.channels[c].bandwidth / 2;
          Channel child1(newBw), child2(newBw);

          vector<Connection> connections = sp.channels[c].connections;

          for (Connection &conn : connections) {
            if (dist(whatever) % 2) {  // TODO: perigoso aqui
              child1.connections.emplace_back(conn);
            } else {
              child2.connections.emplace_back(conn);
            }
          }

          sp.channels[c] = child1;
          sp.channels.insert(sp.channels.begin() + c,
                             child2);  // TODO: it was .end(). Did it right?
        }
      }
    }
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
            spectrum.channels[curr].bandwidth ==
                spectrum.channels[curr + 1].bandwidth) {
          Channel merged = spectrum.channels[curr];
          merged.bandwidth *= 2;

          for (const auto &conn : spectrum.channels[curr + 1].connections) {
            merged.connections.emplace_back(conn);
          }

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

int random_number(int max_num = MAX_CONN) { return dist(whatever) % max_num; }

void K_AddDrop(Solution &sol, int K) {
  K = min(K, int(sol.slots[0].spectrums[3].channels[0].connections.size()));
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
    int a = random_number(sol(t).spectrums.size());
    int b = random_number(sol(t, a).channels.size());

    if (sol(t, a, b).connections.empty()) continue;

    Channel &ch = sol(t, a, b);
    int z = random_number(ch.connections.size());
    Connection conn = ch.connections[z];

    reinsert(sol, conn, make_tuple(t, a, b), zeroChannel, true);
    k++;
  }

  K_AddDrop(sol, K);
}

void perturbation3(Solution &inc, Solution &sol20, int K) {
  vector<pair<Connection, tuple<int, int, int>>> conns;
  for (int t = 0; t < inc.slots.size(); t++)
    for (int s = 0; s < inc(t).spectrums.size(); s++)
      for (int c = 0; c < inc(t, s).channels.size(); c++)
        for (int l = 0; l < inc(t, s, c).connections.size(); l++)
          conns.emplace_back(
              make_pair(inc(t, s, c).connections[l], make_tuple(-1, -1, -1)));

  for (int t = 0; t < sol20.slots.size(); t++)
    for (int s = 0; s < sol20(t).spectrums.size(); s++)
      for (int c = 0; c < sol20(t, s).channels.size(); c++)
        for (int l = 0; l < sol20(t, s, c).connections.size(); l++) {
          for (int x = 0; x < conns.size(); x++) {
            if (conns[x].first.id == sol20(t, s, c).connections[l].id)
              conns[x].second = make_tuple(t, s, c);
          }
        }

  sort(conns.rbegin(), conns.rend(),
       [](const pair<Connection, tuple<int, int, int>> &a,
          const pair<Connection, tuple<int, int, int>> &b) {
         return definitelyLessThan(a.first.violation, b.first.violation);
       });

  K = min(K, int(conns.size()));
  for (int k = 0; k < K; k++) {
    int t = get<0>(conns[k].second);
    int a = get<1>(conns[k].second);
    int b = get<2>(conns[k].second);
    reinsert(sol20, conns[k].first, make_tuple(t, a, b), zeroChannel, true);
  }
  K_AddDrop(sol20, K);
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

  for (int t = 0; t < sol.slots.size(); t++) {
    for (int s = 0; s < sol(t).spectrums.size(); s++) {
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        chanOF[t][s][c] = sol(t, s, c).throughput;
      }
    }
  }
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

int compareObjectives(const int lhs, const int rhs) {
  if (definitelyGreaterThan(lhs, rhs)) return -1;
  return essentiallyEqual(lhs, rhs) ? 0 : 1;
}

int compareObjectives(const Solution &lhs, const Solution &rhs) {
  return compareObjectives(lhs.throughput, rhs.throughput);
}

inline void removeAllOccurrences(Solution &sol, int id) {
  for (int t = 0; t < sol.slots.size(); t++) {
    for (int s = 0; s < sol(t).spectrums.size(); s++) {
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
  }
}

inline void addEverywhere(Solution &sol, int id) {
  for (int t = 0; t < sol.slots.size(); t++) {
    for (int s = 0; s < sol(t).spectrums.size(); s++) {
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        Channel &ch = sol(t, s, c);
        ch = insertInChannel(ch, id);
      }
    }
  }
}

Solution local_search(Solution &mult, Solution &sol20) {
  bool improved = false;
  auto dp = [](const Solution &s) -> double { return s.throughput; };

  do {
    improved = false;

    vector<int> target;
    for (int i = 0; i < n_connections; i++) target.emplace_back(i);

    for (int i : target) {
      Solution mult_clean(mult);
      removeAllOccurrences(mult_clean, i);

      Solution mult_cont(mult_clean);
      addEverywhere(mult_cont, i);

      double bestOF = -1.0;
      ti3 best_ch = make_tuple(-1, -1, -1);

      setDP(mult_clean);
      double cleanOF = calcDP(mult);

      for (int t = 0; t < sol20.slots.size(); t++) {
        for (int s = 0; s < sol20(t).spectrums.size(); s++) {
          for (int c = 0; c < sol20(t, s).channels.size(); c++) {
            setDP(mult_clean, t, s);
            double AA = calcDP(mult, t, s);

            setDP(mult_clean, t, s);
            int c_ch = c;
            while (c_ch != -1) {
              chanOF[t][s][c_ch] = mult_cont(t, s).channels[c_ch].throughput;
              c_ch = parent[t][s][c_ch];
            }

            double BB = calcDP(mult, t, s);
            double result = cleanOF - AA + BB;

            // setDP(mult_clean);
            // int c_ch = c;
            // while (c_ch != -1) {
            //     chanOF[t][s][c_ch] =
            //         mult_cont(t, s).channels[c_ch].throughput;
            //     c_ch = parent[t][s][c_ch];
            // }
            //
            // double result = calcDP(mult);

            if (compareObjectives(result, bestOF) < 0) {
              bestOF = result;
              best_ch = make_tuple(t, s, c);
            }
          }
        }
      }

      int better = compareObjectives(bestOF, dp(mult));
      printf("%lf %lf\n", bestOF, mult.throughput);
      if (better < 0 && best_ch != make_tuple(-1, -1, -1)) {
        improved = true;

        int aux = dp(mult);
        mult = mult_clean;
        mult.throughput = aux;

        int nt = get<0>(best_ch);
        int nsp = get<1>(best_ch);
        int nch = get<2>(best_ch);
        while (nch != -1) {
          mult.slots[nt].spectrums[nsp].channels[nch] =
              insertInChannel(mult.slots[nt].spectrums[nsp].channels[nch], i);
          nch = parent[nt][nsp][nch];
        }
      }

      printf("connection %d better %d\n", i, improved ? 1 : 0);
    }

  } while (improved);

  return reconstruct_sol(mult);
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

vector<Channel> split(Channel toSplit) {
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
  return ret;
}

struct ChanComparer {
  bool operator()(const Channel &a, const Channel &b) const {
    return definitelyLessThan(a.violation, b.violation);
  }
};

void perturbation2(Solution &sol20, int K) {
  vector<Connection> conns;
  for (const auto &s : sol20.slots)
    for (const auto &sp : s.spectrums)
      for (const auto &ch : sp.channels)
        for (const auto &link : ch.connections)
          if (definitelyGreaterThan(link.violation, 0.0))
            conns.emplace_back(link);

  sort(conns.rbegin(), conns.rend(),
       [](const Connection &a, const Connection &b) {
         return definitelyLessThan(a.violation, b.violation);
       });

  K = min(int(conns.size()), K);
  int k = 0;
  while (k < K) {
    Connection to_rmv = conns[0];

    bool ok = false;
    int t, s, c, l;

    for (t = 0; t < sol20.slots.size(); t++) {
      for (s = 0; s < sol20(t).spectrums.size(); s++) {
        for (c = 0; c < sol20(t, s).channels.size(); c++) {
          for (l = 0; l < sol20(t, s, c).connections.size(); l++) {
            if (sol20(t, s, c).connections[l].id == to_rmv.id) {
              ok = true;
              goto fora;
            }
          }
        }
      }
    }

  fora:
    assert(ok);
    reinsert(sol20, to_rmv, make_tuple(t, s, c), zeroChannel, true);

    conns.erase(conns.begin());
    ++k;
  }
  K_AddDrop(sol20, K);
}

inline void perturbation(Solution &sol, int kkmul) {
  int rnd = random_number(2);
  if (rnd)
    K_AddDrop(sol, kkmul);
  else
    K_RemoveAndInserts(sol, kkmul);

  computeThroughput(sol);
  fix_channels(sol);
}

bool checkVio(Solution &sol) {
  double m_ = numeric_limits<double>::min();

  for (int t = 0; t < sol.slots.size(); t++) {
    for (int s = 0; s < sol(t).spectrums.size(); s++) {
      for (int c = 0; c < sol(t, s).channels.size(); c++) {
        for (const auto &conn : sol(t, s, c).connections) {
          m_ = max(m_, conn.violation);
        }
      }
    }
  }

  assert(approximatelyEqual(sol.violation, m_));
  return approximatelyEqual(sol.violation, m_);
}

Solution CH_VRBSP() {
  Solution ret(init_conf, 0.0, true);
  vector<int> links;
  for (int i = 0; i < n_connections; i++) links.emplace_back(i);

  shuffle(links.begin(), links.end(), whatever);

  for (int i = 0; i < n_connections; ++i) {
    Solution tmp = ret;

    for (int s = 0; s < ret(0).spectrums.size(); ++s)
      for (int c = 0; c < ret(0, s).channels.size(); ++c) {
        Solution s1, s2, s3;

        vector<Channel> channels = split(ret(0, s, c));

        Solution *aux[] = {&tmp, &s1, &s2, &s3};
        tmp = **max_element(begin(aux), end(aux), [](Solution *a, Solution *b) {
          return definitelyLessThan(a->throughput, b->throughput);
        });
      }

    ret = tmp;
  }

  return ret;
}

#endif