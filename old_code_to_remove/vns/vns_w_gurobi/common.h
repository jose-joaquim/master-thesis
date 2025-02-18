#ifndef COMMON_H
#define COMMON_H

#include <cassert>
#include <chrono>
#include <cstdio>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "gurobi_c++.h"

const int MAX_CONN = 2048;
const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-9;

using namespace std;
using namespace std::chrono;

using ii = pair<int, int>;
using ti3 = tuple<int, int, int>;
using mt3 = map<tuple<int, int, int>, GRBVar>;
using mt2 = map<tuple<int, int>, GRBVar>;

extern int N;
extern high_resolution_clock::time_point start;
extern double maximumTime;

// static const int C = 45;

static vector<double> GMM, BM;

struct Connection {
  int id;
  double throughput;
  double interference;
  double SINR;
  double distanceSR;
  double violation;

  Connection(int id, double throughput, double interference, double distanceSR);
  Connection(int id);

  bool operator<(const Connection &other) const;
  bool operator>(const Connection &other) const;
};

struct Channel {
  double throughput;
  double interference;
  double violation;
  int bandwidth;
  vector<Connection> connections;

  Channel(double throughput, double interference, double violation,
          int bandwidth, const vector<Connection> connections);
  Channel(int bandwidth);
  Channel();

  Connection &operator()(int l);
  bool operator<(const Channel &other) const;
  bool operator>(const Channel &other) const;
};

struct Spectrum {
  int maxFrequency;
  int usedFrequency;
  double interference;
  vector<Channel> channels;

  Spectrum(int maxFrequency, int usedFrequency, const vector<Channel> channels);
  Spectrum();

  Channel &operator()(int c);
  Connection &operator()(int c, int l);
};

struct TimeSlot {
  vector<Spectrum> spectrums;
  double interference;
  double throughput;

  TimeSlot(const vector<Spectrum> sp);
  TimeSlot();

  Spectrum &operator()(int s);
  Channel &operator()(int s, int c);
  Connection &operator()(int s, int c, int l);
};

class Solution {
 public:
  vector<TimeSlot> slots;
  double throughput;
  double violation;

  Solution(const vector<Spectrum> sp, double tot, bool flag);
  Solution(const vector<TimeSlot> &ts);
  Solution();

#ifdef MDVRBSP
  bool operator<(const Solution &o1) const;
  bool operator>(const Solution &o1) const;
#else
  bool operator<(const Solution &o1) const;
  bool operator>(const Solution &o1) const;
#endif

  TimeSlot &operator()(int t);
  Spectrum &operator()(int t, int s);
  Channel &operator()(int t, int s, int c);
  Connection &operator()(int t, int s, int c, int l);
};

void init(void);

Solution vns(Solution, string);

Solution vns(string);

bool is_feasible(const Solution &, bool);

int bwIdx(int bw);

bool approximatelyEqual(double a, double b, double epsilon = EPS);

bool essentiallyEqual(double a, double b, double epsilon = EPS);

bool definitelyGreaterThan(double a, double b, double epsilon = EPS);

bool definitelyLessThan(double a, double b, double epsilon = EPS);

double computeThroughput(Solution &curr, bool force = false);

bool reinsert(Solution &sol, Connection conn, ti3 from, ti3 to,
              bool force = false);

Channel deleteFromChannel(const Channel &channel, int idConn);

Channel insertInChannel(Channel newChannel, int idConn);

double computeConnectionThroughput(Connection &conn, int bandwidth);

void rmvUnusedTS(Solution &sol);

bool stop(void);

double distance(double X_si, double Y_si, double X_ri, double Y_ri);

double gammaToBeta(double gmm, int bw);

double convertDBMToMW(double _value);

void distanceAndInterference();

void convertTableToMW(const vector<vector<double>> &_SINR,
                      vector<vector<double>> &__SINR);

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static vector<Spectrum> init_conf;
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#ifndef MDVRBSP
void read_data();
#else
void read_data(vector<double> &GMA);
#endif

void initTimeSlot();

void print_objective_to_file(const Solution &sol, FILE *file = nullptr);

void print_solution_to_gurobi(const Solution &sol, FILE *file);

bool fix_channels(Solution &sol);

void K_RemoveAndInserts(Solution &sol, int K);

void K_AddDrop(Solution &sol, int K);

double solve(int k, int i, int j);

int compareObjectives(const int lhs, const int rhs);

int compareObjectives(const Solution &lhs, const Solution &rhs);

optional<vector<Channel>> split(Channel toSplit, vector<double> *GMA = nullptr);

void computeChannelsThroughput(vector<Channel> &channels);

Solution perturbation(Solution &sol, int kkmul);

Solution CH_VRBSP();

Solution ls_math_optimizer(const Solution &);

vector<pair<int, int>> run_math_solver(vector<int> &links, int spectrum_idx);

#endif
