#ifndef MASTER_BASIC_H_
#define MASTER_BASIC_H_

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <map>
#include <nlohmann/json.hpp>
#include <random>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

const int C = 45, M = 12, BAND = 4;
const int MAX_CONN = 2048;
const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-6;

using namespace std;
using namespace std::chrono;
using json = nlohmann::json;
using hrc = std::chrono::high_resolution_clock;
using dii = pair<double, pair<int, int>>;
using misi = map<int, unordered_set<int>>;

#if defined(USE_MATH_SOLVER) || defined(USE_VRBSP_IP) || defined(USE_MDVRBSP_IP)
#include "gurobi_c++.h"
#endif

using ii = pair<int, int>;

extern FILE *progress_file;

extern int N, T;
extern double alfa, NOI, powerSender;
extern double receivers[MAX_CONN][2], senders[MAX_CONN][2];
extern double AFF[MAX_CONN][MAX_CONN];
extern double DM[MAX_CONN][MAX_CONN];
extern vector<vector<double>> DR, SINR, B;
extern vector<unordered_set<int>> C_b;

extern vector<double> GMM, BM;

extern random_device rd;
extern default_random_engine whatever;

extern vector<vector<int>> overlap;

struct Connection {
  int id;
  int mcs;
  double throughput;
  double interference;
  double SINR;
  double distanceSR;
  double violation;

  Connection(int, double, double, double);
  Connection(int);

  bool operator<(const Connection &) const;
  bool operator>(const Connection &) const;
};

struct Channel {
  int ix;
  double throughput;
  double interference;
  double violation;
  int bandwidth;
  int parent;
  int child[2];
  bool in_solution;

  vector<Connection> connections;

  Channel(int, const vector<Connection> &);
  Channel();

  Connection &operator()(int);
  bool operator<(const Channel &) const;
  bool operator>(const Channel &) const;
};

struct TimeSlot {
  vector<Channel> channels;
  double interference;
  double throughput;
  double used_bw;
  unordered_map<int, int> ub_channels;

  TimeSlot(const vector<Channel> &);
  TimeSlot();

  Channel &operator()(int);
  Connection &operator()(int, int);
};

class Solution {
public:
  unordered_set<int> links_not_inserted;
  vector<TimeSlot> slots;
  double throughput_;
  double violation;

  Solution(const vector<TimeSlot> &);
  Solution();

  bool insert_connnection(int, int, int);
  int remove_connection(int);

  double optimal_partitioning();

  Solution get_optimal_solution();
  void get_optimal_solution(Solution &, int, int, bool);

  vector<Connection> get_scheduled_connections() const;
  double decode(std::vector<double>) const;

  bool operator<(const Solution &) const;
  bool operator>(const Solution &) const;
  bool operator==(const Solution &) const;

  TimeSlot &operator()(int t);
  Channel &operator()(int t, int c);
  Connection &operator()(int t, int c, int l);
};

void to_json(json &, const Solution &);

void init(void);

vector<Solution> vns(string);

bool is_feasible(const Solution &, bool);

int bwIdx(int);

int bToIdx(int);

bool approximatelyEqual(double, double, double = EPS);

bool essentiallyEqual(double, double, double = EPS);

bool definitelyGreaterThan(double, double, double = EPS);

bool definitelyLessThan(double, double, double = EPS);

bool GreaterOrEqual(double, double);

double computeThroughput(Solution &, bool = false);

void reinsert(Solution &, Connection, ii, ii);

optional<Channel> deleteFromChannel(Channel, int);

optional<Channel> insertInChannel(Channel, int, bool = false);

double computeConnectionThroughput(Connection &conn, int bandwidth);

void rmvUnusedTS(Solution &sol);

bool stop(high_resolution_clock::time_point, double);

double distance(double X_si, double Y_si, double X_ri, double Y_ri);

double gammaToBeta(double gmm, int bw);

double convertDBMToMW(double _value);

void distanceAndInterference();

void convertTableToMW(const vector<vector<double>> &_SINR,
                      vector<vector<double>> &__SINR);

void read_data();

void initTimeSlot();

int random_number(int = MAX_CONN);

void print_objective_to_file(const Solution &sol, FILE *file = nullptr);

vector<string> gurobi_sol(const Solution &);

Solution reconstruct_sol(Solution &curr);

Solution multipleRepresentation(Solution ret);

Solution convertTo20MhzSol(Solution &);

bool allChannels20MHz(const Solution &sol);

void K_RemoveAndInserts(Solution &sol, int K);

void K_AddDrop(Solution &sol, int K);

double calcDP(Solution &, int, int, int = 0);

int compareObjectives(const int lhs, const int rhs);

int compareObjectives(const Solution &lhs, const Solution &rhs);

void removeAllOccurrences(Solution &sol, int id);

double optimal_partitioning_global(Solution &, bool = true);

Solution local_search(const Solution &);

optional<vector<Channel>> split(Channel toSplit, vector<double> *GMA = nullptr);

void computeChannelThroughput(Channel &);

Solution perturbation(Solution &sol, int kkmul);

Solution CH_VRBSP();

Solution CH_MDVRBSP();

Solution local_search(const Solution &);

Solution rebuild_solution(Solution);

int necessaryBw(double);

bool canTransmitUsingChannel(int, int);

bool canTransmitUsingBandwidth(int, int, int);

int cToBIdx(int);

int idxToB(int);

int ChannelIdToBandwidth(int);

void fix_channels(Solution &);

double buildVRBSPSolution(vector<double>, vector<int>,
                          optional<string> = std::nullopt);

void label_channels(Solution &);

#if defined(USE_MATH_SOLVER) || defined(USE_VRBSP_IP) || defined(USE_MDVRBSP_IP)
void var_Iij(GRBModel *, GRBVar[][C], misi &);

void var_z(GRBModel *, GRBVar[][C], misi &);

void var_I(GRBModel *, GRBVar[], misi &);

void var_t(GRBModel *, GRBVar *, vector<int> &, misi &);

void var_x(GRBModel *, GRBVar[][C], misi &);

void var_y(GRBModel *, GRBVar[][BAND][M], misi &);

void unique(GRBModel *, GRBVar[][C]);

void couple(GRBModel *, GRBVar[][BAND][M], GRBVar[][C]);

void ch_overlap(GRBModel *, GRBVar[][C], GRBVar[][C]);

void interch(GRBModel *, GRBVar[][C], GRBVar[][C]);

void bigG(GRBModel *, GRBVar[], GRBVar[][C], GRBVar[][C]);

void bigL(GRBModel *, GRBVar[], GRBVar[][C], GRBVar[][C]);

void sinr(GRBModel *, GRBVar[], GRBVar[][BAND][M]);

Solution vrbsp(misi &);
#endif

#endif
