#include "basic.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <type_traits>

#if defined(USE_MATH_SOLVER)
#include "gurobi_c++.h"

using mt3 = map<tuple<int, int, int>, GRBVar>;
#endif

FILE *progress_file;

int N, T, cnt;
double alfa, NOI, powerSender;
double receivers[MAX_CONN][2], senders[MAX_CONN][2];
double AFF[MAX_CONN][MAX_CONN];
double DM[MAX_CONN][MAX_CONN];
vector<vector<double>> DR, SINR, B;

random_device rd;
mt19937 gen{rd()};
uniform_int_distribution<> dist{0, MAX_CONN};
default_random_engine whatever = default_random_engine{rd()};
vector<double> GMM, BM;

vector<vector<int>> overlap = {
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

Connection::Connection(int id, double throughput, double interference,
                       double distanceSR)
    : id(id), throughput(throughput), interference(interference),
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

Channel::Channel(int bandwidth, const vector<Connection> &conns)
    : bandwidth(bandwidth), connections(conns) {
  throughput = 0.0;
  interference = 0.0;
  violation = 0.0;
  parent = -1;
  in_solution = false;
  memset(child, -1, sizeof child);
}

Channel::Channel() : Channel(0.0, {}) {}

Connection &Channel::operator()(int l) { return connections[l]; }

TimeSlot::TimeSlot() {
  interference = 0.0;
  throughput = 0.0;
  ub_channels = {{160, 0}, {80, 0}, {40, 0}, {20, 0}};
}

TimeSlot::TimeSlot(const vector<Channel> &chs) : channels(chs) {
  interference = 0.0;
  throughput = 0.0;

  ub_channels = {{160, 0}, {80, 0}, {40, 0}, {20, 0}};

  for (const Channel &ch : channels)
    ub_channels[ch.bandwidth] += 1;
}

Channel &TimeSlot::operator()(int c) { return channels[c]; }

Connection &TimeSlot::operator()(int c, int l) {
  return channels[c].connections[l];
}

vector<Connection> Solution::get_scheduled_connections() const {
  vector<Connection> ret;
  for (const TimeSlot &ts : slots)
    for (const Channel &ch : ts.channels)
      for (const Connection &conn : ch.connections)
        ret.push_back(conn);

  return ret;
}

void to_json(json &json, const Solution &sol) {
  vector<pair<int, vector<int>>> channels;

  for (const Channel &ch : sol.slots[0].channels) {
    vector<int> conns;
    for (const Connection &conn : ch.connections)
      if (definitelyGreaterThan(conn.throughput, 0.0))
        conns.push_back(conn.id);

    channels.push_back({ch.bandwidth, conns});
  }

  json["throughput"] = sol.throughput_;
  json["channels"] = channels;
}

Solution::Solution(const vector<TimeSlot> &ts) : slots(ts) {
  throughput_ = 0.0;
  violation = 0.0;
}

Solution::Solution() {
  slots = vector<TimeSlot>();
  throughput_ = 0.0;
  violation = 0.0;
}

#ifdef MDVRBSP
// TODO: update here
bool Solution::operator<(const Solution &o1) const {
  return violation < o1.violation;
}

bool Solution::operator>(const Solution &o1) const {
  return !Solution::operator<(o1) && !Solution::operator==(o1);
}
#else
bool Solution::operator<(const Solution &o1) const {
  double diff = throughput_ - o1.throughput_;
  return definitelyLessThan(diff, 0.0) && !approximatelyEqual(diff, 0.0);
}

bool Solution::operator>(const Solution &o1) const {
  return !Solution::operator<(o1) && !Solution::operator==(o1);
}
#endif

bool Solution::operator==(const Solution &o1) const {
  return approximatelyEqual(throughput_, o1.throughput_);
}

TimeSlot &Solution::operator()(int t) { return slots[t]; }

Channel &Solution::operator()(int t, int c) { return slots[t].channels[c]; }

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
  return hypot((X_si - X_ri), (Y_si - Y_ri));
}

void fix_channels(Solution &sol) {
  vector<int> infeas_conns;
  for (const Channel &ch : sol(0).channels)
    for (const Connection &conn : ch.connections)
      if (approximatelyEqual(0.0, conn.throughput))
        infeas_conns.push_back(conn.id);

  // if (!infeas_conns.empty())
  //   printf("removing %lu\n", infeas_conns.size());

  for (const int ix : infeas_conns)
    sol.remove_connection(ix);

  computeThroughput(sol);
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

bool SINRfeas(const Connection &con, int bw) {
  return GreaterOrEqual(con.SINR, B[con.id][bToIdx(bw)]);
}

double computeConnectionThroughput(Connection &i, int bw) {
  if (bw == 0)
    return 0.0;

  if (approximatelyEqual(i.interference, 0.0, EPS)) {
    i.SINR = SINR[11][bToIdx(bw)];
    i.throughput = DR[11][bToIdx(bw)];
  } else {
    i.SINR = AFF[i.id][i.id] / (i.interference + NOI);

    int mcs = 11;
    while (mcs >= 0 && definitelyLessThan(i.SINR, SINR[mcs][bToIdx(bw)]))
      mcs--;

    if (mcs < 0) {
      i.throughput = 0.0;
    } else
      i.throughput = DR[mcs][bToIdx(bw)];
  }

  return i.throughput;
}

optional<Channel> deleteFromChannel(Channel channel, int idConn) {
  for (int j = 0; const Connection &conn : channel.connections)
    if (conn.id == idConn) {
      swap(channel.connections[j], channel.connections.back());
      channel.connections.pop_back();
      break;
    } else {
      j++;
    }

  channel.throughput = 0.0;
  channel.violation = 0.0;
  for (Connection &conn : channel.connections) {
    conn.interference -= AFF[idConn][conn.id];

    computeConnectionThroughput(conn, channel.bandwidth);
    channel.throughput += conn.throughput;
  }

  return channel;
}

optional<Channel> insertInChannel(Channel new_c, int i_idx, bool force) {
  Connection new_conn(i_idx, 0.0, 0.0, DM[i_idx][i_idx]);

  for (auto &j : new_c.connections) {
    j.interference += AFF[i_idx][j.id];
    new_conn.interference += AFF[j.id][i_idx];
  }

  new_c.throughput = 0.0;
  new_c.connections.emplace_back(new_conn);
  for (auto &connection : new_c.connections) {
    double ans = computeConnectionThroughput(connection, new_c.bandwidth);
    if (!force && !SINRfeas(connection, new_c.bandwidth))
      return nullopt;

    new_c.throughput += connection.throughput;
  }

  return new_c;
}

optional<vector<Channel>> split(Channel toSplit, vector<double> *GMA) {
  vector<Channel> ret(2, Channel(toSplit.bandwidth / 2, {}));

  vector<Connection> reinsert = toSplit.connections;
  sort(reinsert.rbegin(), reinsert.rend());

  for (int i = 0; i < reinsert.size(); i++) {
    optional<Channel> bestChannel;
    int best_c = -1;

    for (int c = 0; c < ret.size(); c++) {
      optional<Channel> inserted = insertInChannel(ret[c], reinsert[i].id);

      if (inserted && inserted->throughput > bestChannel->throughput) {
        bestChannel = inserted;
        best_c = c;
      }
    }

    if (!bestChannel)
      return nullopt;

    assert(best_c >= 0);
    ret[best_c] = *bestChannel;
  }

  if (GMA != nullptr)
    for (const auto &ch : ret)
      for (const auto &conn : ch.connections)
        if (definitelyLessThan(conn.throughput, (*GMA)[conn.id]))
          return nullopt;

  return ret;
}

optional<Channel> try_insert(int conn_id, optional<Channel> ch) {
  if (definitelyGreaterThan(GMM[conn_id], DR[11][bToIdx(ch->bandwidth)]))
    return nullopt;

  ch = insertInChannel(*ch, conn_id);
  if (!ch)
    return nullopt;

  for (Connection &conn : ch->connections)
    if (definitelyLessThan(conn.SINR, B[conn.id][bToIdx(ch->bandwidth)]))
      return nullopt;

  return ch;
}

bool can_split(const Channel &ch) {
  if (ch.bandwidth <= 20)
    return false;

  for (const auto &conn : ch.connections)
    if (definitelyLessThan(DR[11][bToIdx(ch.bandwidth / 2)], GMM[conn.id]))
      return false;

  return true;
}

#if defined(USE_DECODER)
int convertChromoBandwidth(double value) {
  if (value < .25)
    return 20;
  else if (value < .5)
    return 40;
  else if (value < .75)
    return 80;

  return 160;
}

unordered_map<int, int> ub_channels = {{160, 2}, {80, 6}, {40, 12}, {20, 25}};

int insertFreeChannel(Solution &sol, int conn, int band,
                      vector<double> &variables) {
  int alternative_bw = -1;
  int largest_free_bw = 0;
  for (TimeSlot &ts : sol.slots) {
    int free_bw = 500 - ts.used_bw;
    if (free_bw >= band && sol(0).ub_channels[band] + 1 < ub_channels[band]) {
      Channel newChannel = Channel(band, {});
      newChannel = *insertInChannel(newChannel, conn, true);
      ts.channels.emplace_back(newChannel);
      ts.used_bw += band;

      sol(0).ub_channels[band] += 1;
      sol.throughput_ += newChannel.throughput;
      return band;
    } else if (free_bw >= largest_free_bw) {
      largest_free_bw = free_bw;

      if (largest_free_bw >= 80 &&
          sol(0).ub_channels[band] + 1 < ub_channels[80])
        alternative_bw = 80;
      else if (largest_free_bw >= 40 &&
               sol(0).ub_channels[band] + 1 < ub_channels[40])
        alternative_bw = 40;
      else if (largest_free_bw >= 20 &&
               sol(0).ub_channels[band] + 1 < ub_channels[20])
        alternative_bw = 20;
    }
  }

  assert(alternative_bw > 0);
  sol(0).channels.emplace_back(Channel(alternative_bw, {conn}));
  sol(0).used_bw += alternative_bw;
  sol(0).ub_channels[alternative_bw] += 1;
  computeThroughput(sol);

  if (alternative_bw != band) {
    switch (alternative_bw) {
    case 20:
      variables[(conn * 2) + 1] = (0 + 0.25) / 2.0;
      break;
    case 40:
      variables[(conn * 2) + 1] = (0.5 + 0.25) / 2.0;
      break;
    case 80:
      variables[(conn * 2) + 1] = (0.75 + 0.5) / 2.0;
      break;
    case 160:
      variables[(conn * 2) + 1] = (1.0 + 0.75) / 2.0;
      break;
    default:
      exit(77);
    }
  }

  return alternative_bw;
}

void insertBestChannel(Solution &sol, int conn, int band,
                       vector<double> &variables) {
  pair<int, int> nCh = {-1, -1};
  double max_throughput = sol.throughput_;
  Channel seilaaaa;
  for (int t = 0; t < int(sol.slots.size()); t++) {
    for (int c = 0; c < sol(t).channels.size(); c++) {
      Channel tmp_channel = *insertInChannel(sol(t, c), conn, true);
      double update_throughput =
          sol.throughput_ - sol(t, c).throughput + tmp_channel.throughput;

      if (definitelyGreaterThan(update_throughput, max_throughput)) {
        max_throughput = update_throughput;
        seilaaaa = tmp_channel;
        nCh = {t, c};
      }
    }
  }

  if (nCh == make_pair(-1, -1))
    return;

  const auto &[ts, ch] = nCh;
  Channel old_channel = sol(ts, ch);
  sol(ts, ch) = seilaaaa;
  sol.throughput_ = max_throughput;

  if (seilaaaa.bandwidth != band) {
    switch (seilaaaa.bandwidth) {
    case 20:
      variables[(conn * 2) + 1] = (0 + 0.25) / 2.0;
      break;
    case 40:
      variables[(conn * 2) + 1] = (0.5 + 0.25) / 2.0;
      break;
    case 80:
      variables[(conn * 2) + 1] = (0.75 + 0.5) / 2.0;
      break;
    case 160:
      variables[(conn * 2) + 1] = (1.0 + 0.75) / 2.0;
      break;
    default:
      exit(77);
    }
  }
}

double buildVRBSPSolution(vector<double> variables, vector<int> permutation,
                          optional<string> out_file) {
  // printf("start\n");
  int totalSpectrum = 500, used_bw = 0;

  // Channel c160(160, {}), c80(80, {}), c40(40, {}), c20(20, {});
  TimeSlot init(vector<Channel>{});
  init.used_bw = 0.0;
  Solution sol(vector<TimeSlot>{init});

  // First, insert in free channels
  int idx = 0;
  while (idx < permutation.size() && used_bw < totalSpectrum) {
    int connection = permutation[idx] / 2;
    int bandWidth = convertChromoBandwidth(variables[permutation[idx] + 1]);
    used_bw += insertFreeChannel(sol, connection, bandWidth, variables);
    idx += 1;
  }

  computeThroughput(sol);

  // Second, insert in the best channels
  while (idx < permutation.size()) {
    int connection = permutation[idx] / 2;
    int bandWidth = convertChromoBandwidth(variables[permutation[idx] + 1]);
    insertBestChannel(sol, connection, bandWidth, variables);
    idx += 1;
  }

  fix_channels(sol);

  Solution seila = convertTo20MhzSol(sol);
  Solution best_multiple = multipleRepresentation(seila);
  double dp_of = optimal_partitioning_global(best_multiple, false);
  if (!(definitelyGreaterThan(dp_of, sol.throughput_) ||
        approximatelyEqual(dp_of, sol.throughput_)))
    printf("%lf %lf\n", dp_of, sol.throughput_);

  int cnt_160 = 0, cnt_80 = 0, cnt_40 = 0, cnt_20 = 0;

  for (const Channel &ch : sol(0).channels)
    if (ch.bandwidth == 160)
      cnt_160++;
    else if (ch.bandwidth == 80)
      cnt_80++;
    else if (ch.bandwidth == 40)
      cnt_40++;
    else
      cnt_20++;

  assert(cnt_160 <= 2);
  assert(cnt_80 <= 6);
  assert(cnt_40 <= 12);
  assert(cnt_20 <= 25);

  if (out_file.has_value()) {
    std::ofstream o(*out_file);
    json seila = json(sol);
    o << seila.dump(4);
    o.close();
  }

  return -1.0 * sol.throughput_;
}

double Solution::decode(std::vector<double> variables) const {
  vector<pair<double, int>> ranking;
  for (int i = 0; i < variables.size(); i += 2)
    ranking.emplace_back(variables[i], i);

  sort(ranking.begin(), ranking.end());

  vector<int> permutation;
  for (int i = 0; i < ranking.size(); i++)
    permutation.emplace_back(ranking[i].second);

  double fitness = buildVRBSPSolution(variables, permutation);
  return fitness;
}
#endif

#if defined(USE_MDVRBSP_IP)
Solution CH_MDVRBSP() {
  // TODO: need to fix because I included timeslots
  Channel ch160(160);
  Channel ch80(80);
  Channel ch20(20);
  vector<Channel> chs{ch160, ch160, ch80, ch20};
  TimeSlot dummy_ts(chs);
  Solution ret({dummy_ts});

  vector<int> links(N);
  iota(links.begin(), links.end(), 0);
  shuffle(links.begin(), links.end(), whatever);

  for (int conn : links) {
    for (int t = 0; t < ret.slots.size(); t++) {
      if (conn < t)
        continue;
      for (int c = 0; c < ret(t).channels.size(); c++) {
        auto rst = try_insert(conn, ret(t, c));

        if (rst) {
          ret(t, c) = *rst;

          goto NEXT_CONN;
        }

        if (can_split(ret(t, c))) {
          auto ch_split = split(ret(t, c), &GMM);

          if (!ch_split)
            continue;

          for (int i = 0; i < ch_split->size(); ++i) {
            auto cp_channel = try_insert(conn, (*ch_split)[i]);

            if (cp_channel) {
              swap(ret(t, c), ret(t).channels.back());
              ret(t).channels.pop_back();
              ret(t).channels.emplace_back(*cp_channel);
              ret(t).channels.emplace_back((*ch_split)[i == 0]);
              goto NEXT_CONN;
            }
          }
        }
      }
    }

    {
      TimeSlot new_ts{dummy_ts};
      int bw = new_ts(0).bandwidth;
      assert(bw == 160);
      Connection new_conn{conn};

      new_conn.SINR = 1000000007;
      new_conn.throughput = DR[11][bToIdx(bw)];
      new_ts(0).connections.emplace_back(new_conn);
      ret.slots.emplace_back(new_ts);
      goto NEXT_CONN;
    }

  NEXT_CONN:;
  }

  for (int t = 0; const TimeSlot &ts : ret.slots)
    for (const Channel &ch : ts.channels)
      for (const Connection &conn : ch.connections)
        assert(SINRfeas(conn, ch.bandwidth));

  return ret;
}
#endif

Channel erase_connection_from_channel(Channel ch, int id) {
  Channel ch_new = ch;
  ch_new.connections.clear();
  ch_new.throughput = 0.0;
  for (const Connection &i : ch.connections)
    if (i.id != id)
      ch_new = *insertInChannel(ch_new, i.id, true);

  return ch_new;
}

int Solution::remove_connection(int i) {
  int found = -1;

  for (int ts = 0; ts < slots.size(); ++ts)
    for (int c = 0; c < slots[ts].channels.size(); ++c)
      for (const Connection &l : slots[ts].channels[c].connections)
        if (l.id == i) {
          (*this)(ts, c) = erase_connection_from_channel((*this)(ts, c), i);
          found = c;
        }

  return found;
}

void read_data() {
  int n_spectrums;
  scanf("%d %lf %lf %lf %d", &N, &alfa, &NOI, &powerSender, &n_spectrums);

  for (int i = 0; i < n_spectrums; i++) {
    int s;
    scanf("%d", &s);
  }

  if (NOI != 0) {
    NOI = convertDBMToMW(NOI);
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

  GMM.assign(N, 0.0);
  for (int i = 0; i < N; i++) {
    double bt;
    scanf("%lf", &bt);
#if defined(USE_MDVRBSP_IP) || defined(USE_MDVRBSP_COL_GEN)
    GMM[i] = bt;
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

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double value = powerSender / pow(DM[i][j], alfa);
      AFF[i][j] = value;
    }
  }

  // TODO: find a better value for this big-M
  BM.assign(N, 100.0);
  B.assign(N, vector<double>(4, 0));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; ++j) {
      if (i != j)
        BM[i] += AFF[j][i];

      if (i != j) {
        int bwj = necessaryBw(GMM[j]);
        int bwi = necessaryBw(GMM[i]);

        if (bwj >= bwi)
          BM[i] += AFF[j][i] + 0.0005;
      }
    }

    for (int j = 0; j < 4; j++) {
      double value = gammaToBeta(GMM[i], j);
      assert(definitelyGreaterThan(value, 0.0));
      B[i][j] = value;
    }
  }
}

double gammaToBeta(double gmm, int bw) {
  int mcs = 0;

  while (mcs < 12 && definitelyGreaterThan(gmm, DR[mcs][bw]))
    mcs++;

  if (mcs < 12)
    return SINR[mcs][bw];

  return SINR[11][3] + 1.0;
}

double convertDBMToMW(double _value) {
  double result = 0.0, b;

  b = _value / 10.0;     // dBm dividido por 10
  result = pow(10.0, b); // Converte de DBm para mW

  return result;
}

void convertTableToMW(const vector<vector<double>> &_SINR,
                      vector<vector<double>> &__SINR) {
  double result, b;
  for (int i = 0; i < __SINR.size(); i++) {
    for (int j = 0; j < __SINR[i].size(); j++) {
      if (_SINR[i][j] != 0) {
        b = _SINR[i][j] / 10.0; // dBm divided by 10
        result = pow(10.0, b);  // Convert DBM to mW

        __SINR[i][j] = result;
      } else {
        __SINR[i][j] = 0;
      }
    }
  }
}

double computeThroughput(Solution &curr, bool force) {
  curr.throughput_ = 0.0;

  for (TimeSlot &ts : curr.slots) {
    for (Channel &ch : ts.channels) {
      if (ch.bandwidth == 0)
        continue;

      ch.throughput = 0.0;
      for (Connection &u : ch.connections) {
        u.interference = 0.0;
        u.throughput = 0.0;
        for (Connection &v : ch.connections)
          u.interference += u.id == v.id ? 0.0 : AFF[v.id][u.id];

        ch.throughput += computeConnectionThroughput(u, ch.bandwidth);
      }

      curr.throughput_ += ch.throughput;
    }
  }

  return curr.throughput_;
}

void reinsert(Solution &sol, Connection conn, ii from, ii to) {
  if (from == to)
    return;

  const auto &[ft, fc] = from;
  const auto &[tt, tc] = to;

  if (from != make_pair(-1, -1)) {
    for (int j = 0; const Connection &connection : sol(ft, fc).connections)
      if (conn.id == connection.id) {
        swap(sol(ft, fc).connections[j], sol(ft, fc).connections.back());
        sol(ft, fc).connections.pop_back();
        break;
      } else {
        j++;
      }
  }

  if (to != make_pair(-1, -1))
    sol(tt, tc).connections.push_back(conn);

  // double obj_part_remove = 0;
  // double obj_part_add = 0;
  // if (from != make_pair(-1, -1)) {
  //   Channel &old_chan = sol(ft, fc);
  //   obj_part_remove = old_chan.throughput;
  //
  //   old_chan = *deleteFromChannel(old_chan, conn.id);
  //   obj_part_add = old_chan.throughput;
  // }
  //
  // if (to != make_pair(-1, -1)) {
  //   Channel &new_chan = sol(tt, tc);
  //   obj_part_remove += new_chan.throughput;
  //
  //   new_chan = *insertInChannel(new_chan, conn.id, true);
  //   obj_part_add += new_chan.throughput;
  // }
  //
  // sol.throughput_ = sol.throughput_ + obj_part_add - obj_part_remove;
}

void K_RemoveAndInserts(Solution &sol, int K) {
  int k = 0;

  set<int> set_unscheduled;
  for (int i = 0; i < N; ++i)
    set_unscheduled.insert(i);

  K = min(K, int(sol.get_scheduled_connections().size()));
  // assert(!sol.get_scheduled_connections().empty());
  // for (const Connection &conn : sol.get_scheduled_connections()) {
  //   assert(set_unscheduled.contains(conn.id));
  //   set_unscheduled.erase(conn.id);
  // }

  while (k < K) {
    int t = random_number(sol.slots.size());
    int a = random_number(sol(t).channels.size());

    if (sol(t, a).connections.empty())
      continue;

    Channel &ch = sol(t, a);
    int z = random_number(ch.connections.size());
    Connection conn = ch.connections[z];

    reinsert(sol, conn, {t, a}, {-1, -1});
    // printf("removed schedule %d %d %d\n", t, a, conn.id);
    k++;
  }

  K_AddDrop(sol, K);
}

void K_AddDrop(Solution &sol, int K) {
  set<int> set_unscheduled;
  for (int i = 0; i < N; ++i)
    set_unscheduled.insert(i);

  for (const TimeSlot &st : sol.slots)
    for (const Channel &ch : st.channels)
      for (const Connection &conn : ch.connections) {
        assert(set_unscheduled.contains(conn.id));
        set_unscheduled.erase(conn.id);
      }

  vector<int> unscheduled(set_unscheduled.begin(), set_unscheduled.end());
  K = min(K, int(unscheduled.size()));
  for (int i = 0; i < K; ++i) {
    int ix = random_number(unscheduled.size());
    int conn = unscheduled[ix];

    int t = random_number(sol.slots.size());
    int a = random_number(sol(t).channels.size());
    pair<int, int> channelTo = {t, a};

    reinsert(sol, conn, {-1, -1}, channelTo);
    swap(unscheduled[ix], unscheduled.back());
    unscheduled.pop_back();
  }
}

Solution perturbation(Solution &sol, int kkmul) {
  int rnd = dist(gen) % 2;
  if (rnd)
    K_AddDrop(sol, kkmul);
  else
    K_RemoveAndInserts(sol, kkmul);

  computeThroughput(sol);
  fix_channels(sol);

  return sol;
}

Solution convertTo20MhzSol(Solution exps) {
  auto cmp_channels = [](Channel &a, Channel &b) {
    return a.bandwidth < b.bandwidth;
  };
  sort(exps(0).channels.rbegin(), exps(0).channels.rend(), cmp_channels);
  std::uniform_int_distribution<int> distribution(0, 1);
  while (true) {
    Channel *to_replace = 0;
    for (Channel &c : exps(0).channels)
      if (c.bandwidth > 20) {
        to_replace = &c;
        break;
      }

    if (!to_replace)
      break;

    Channel c1(to_replace->bandwidth / 2, {});
    Channel c2 = c1;

    for (const Connection &conn : to_replace->connections) {
      if (distribution(whatever))
        c1 = *insertInChannel(c1, conn.id, true);
      else
        c2 = *insertInChannel(c2, conn.id, true);
    }

    *to_replace = c1;
    exps(0).channels.push_back(c2);
  }

  // Labeling the channels
  for (int t = 0; t < T; ++t)
    for (int c = 0; c < exps.slots[t].channels.size(); ++c)
      exps.slots[t].channels[c].ix = c;

  return exps;
}

int random_number(int max_num) { return dist(gen) % max_num; }

bool stop(hrc::time_point from, double ub) {
  return duration_cast<duration<double>>(hrc::now() - from).count() >= ub;
}

#if defined(USE_VNS_PURE) || defined(USE_VNS_MATH_SOLVER)
Solution multipleRepresentation(Solution ret) {
  for (int ts_idx = 0; TimeSlot & ts : ret.slots) {
    vector<Channel> &chs = ts.channels;
    for (int ix = 0; Channel & ch : ts.channels) {
      ch.in_solution = false;
      ch.parent = -1;
      // ch.ix = ix;
      memset(ch.child, -1, sizeof ch.child);
    }

    int i = 0;
    while (i < chs.size() && chs.size() < 45) {
      if (i + 1 < chs.size() && chs[i].bandwidth < 160 &&
          chs[i].bandwidth == chs[i + 1].bandwidth) {
        Channel merged(2 * chs[i].bandwidth, chs[i].connections);
        for (const Connection &conn : chs[i + 1].connections)
          merged.connections.push_back(conn);

        int p = chs.size();
        chs[i].parent = p;
        chs[i + 1].parent = p;
        merged.child[0] = i;
        merged.child[1] = i + 1;
        merged.ix = p;

        chs.push_back(merged);
        i += 2;
      } else
        i += 1;
    }

    ts_idx++;
  }

  computeThroughput(ret);

  return ret;
}

optional<vector<Channel>> split2(Channel toSplit,
                                 vector<double> *GMA = nullptr) {
  if (toSplit.bandwidth <= 20)
    return nullopt;

  int new_bw = toSplit.bandwidth / 2;
  vector<Connection> conns = toSplit.connections;
  vector<Channel> ret(2, Channel(new_bw, {}));

  for (const auto &con : conns) {
    Channel best(new_bw, {});
    int it_best = -1;

    for (int c = 0; c < ret.size(); ++c) {
      Channel ans = *insertInChannel(ret[c], con.id, true);

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
    return definitelyLessThan(a->throughput_, b->throughput_);
  };

  Channel c160(160, {}), c80(80, {}), c40(40, {}), c20(20, {});
  TimeSlot init(vector<Channel>{c160, c160, c80, c80, c20});
  Solution ret(vector<TimeSlot>{init});
  vector<int> links;
  for (int i = 0; i < N; i++)
    links.emplace_back(i);

  shuffle(links.begin(), links.end(), whatever);

  for (int i = 0; i < N; ++i) {
    Solution tmp = ret;

    for (int c = 0; c < ret(0).channels.size(); ++c) {
      vector<Solution *> candidates{&tmp};
      optional<Channel> c1 = *insertInChannel(ret(0, c), links[i], true);

      Solution s1 = ret, s2 = ret, s3 = ret;
      if (c1.has_value()) {
        s1(0, c) = *c1;
        candidates.push_back(&s1);
      }

      optional<vector<Channel>> opt_channels = split2(ret(0, c));

      if (opt_channels.has_value()) {
        vector<Channel> &channels = *opt_channels;
        Channel c2 = *insertInChannel(channels[0], links[i], true);
        Channel c3 = *insertInChannel(channels[1], links[i], true);

        s2(0, c) = c2;
        s2(0).channels.emplace_back(channels[1]);

        s3(0, c) = c3;
        s3(0).channels.emplace_back(channels[0]);

        candidates.insert(candidates.end(), {&s2, &s3});
      }

      for (Solution *const candidate : candidates)
        computeThroughput(*candidate);

      if (!candidates.empty())
        tmp = **max_element(begin(candidates), end(candidates), compare_sol);
    }

    ret = tmp;
  }

  T = 1;
  return ret;
}

double calcDP(Solution &sol, int t, int c, int depth) {
  // printf("node %d %d %d\n", t, c, depth);
  sol(t, c).in_solution = true;
  double ret = sol(t, c).throughput;

  int child1 = sol(t, c).child[0];
  int child2 = sol(t, c).child[1];
  if (child1 != -1 && child2 != -1) {
    double a = calcDP(sol, t, child1, depth + 1);
    double b = calcDP(sol, t, child2, depth + 1);
    double candidate = a + b;

    if (definitelyGreaterThan(candidate, ret)) {
      ret = candidate;
      sol(t, c).in_solution = false;
    }
  }

  return ret;
}

Solution rebuild_solution(Solution to_rebuild) { return to_rebuild; }

double Solution::optimal_partitioning() {
  // 1. Create multiple representation
  Solution multiple = multipleRepresentation(*this);

  throughput_ = 0.0;
  // 2. Run DP for all maximal channels
  for (int t = 0; t < multiple.slots.size(); ++t)
    for (int c = 0; c < multiple.slots[t].channels.size(); ++c)
      if (multiple(t, c).parent == -1) {
        double ret = calcDP(multiple, t, c);
        throughput_ += ret;
      }

  return throughput_;
}

bool Solution::insert_connnection(int i, int t, int c) {
  optional<Channel> new_channel = *insertInChannel((*this)(t, c), i, true);

  if (!new_channel)
    return false;

  (*this)(t, c) = *new_channel;

  return true;
}

void Solution::get_optimal_solution(Solution &sol, int t, int c, bool remove) {
  if (remove)
    sol(t, c).in_solution = false;

  if (sol(t, c).in_solution)
    remove = true;

  int child1 = sol(t, c).child[0];
  int child2 = sol(t, c).child[1];
  if (child1 != -1 && child2 != -1) {
    get_optimal_solution(sol, t, child1, remove);
    get_optimal_solution(sol, t, child2, remove);
  }
}

Solution Solution::get_optimal_solution() {
  Solution ret;
  double of = 0.0;

  for (int t = 0; t < slots.size(); ++t)
    for (int c = 0; c < slots[t].channels.size(); ++c)
      if ((*this)(t, c).parent == -1)
        get_optimal_solution(*this, t, c, false);

  for (int t = 0; t < slots.size(); ++t) {
    ret.slots.emplace_back(TimeSlot{});
    double sum_bw = 0.0;

    for (int c = 0; c < slots[t].channels.size(); ++c)
      if ((*this)(t, c).in_solution) {
        ret.slots[t].channels.push_back((*this)(t, c));
        of += (*this)(t, c).throughput;
        sum_bw += (*this)(t, c).bandwidth;
      }

    assert(sum_bw == 500.0);
  }

  ret.throughput_ = of;

  return ret;
}

double optimal_partitioning_global(Solution &multiple, bool compute_of) {
  if (compute_of)
    computeThroughput(multiple);

  multiple.throughput_ = 0.0;

  for (int t = 0; t < multiple.slots.size(); ++t)
    for (int c = 0; c < multiple.slots[t].channels.size(); ++c)
      if (multiple(t, c).parent == -1) {
        double ret = calcDP(multiple, t, c);
        multiple.throughput_ += ret;
      }

  return multiple.throughput_;
}

void computeChannelThroughput(Channel &ch) {
  ch.throughput = 0.0;
  for (Connection &u : ch.connections) {
    u.interference = 0.0;
    u.throughput = 0.0;
    for (Connection &v : ch.connections)
      u.interference += u.id == v.id ? 0.0 : AFF[v.id][u.id];

    ch.throughput += computeConnectionThroughput(u, ch.bandwidth);
  }
}

int insert_connection_up_to_root(Solution &sol, int i, int ts, int channel) {
  double highest = -1;
  do {
    sol(ts, channel).connections.push_back(Connection(i));
    computeChannelThroughput(sol(ts, channel));
    channel = sol(ts, channel).parent;

    if (channel != -1)
      highest = channel;

  } while (channel != -1);

  return highest;
}

void erase_connection_up_to_root(Solution &sol, int i, int ts, int channel) {
  do {
    bool deletion_ok = false;
    vector<Connection> &conns = sol(ts, channel).connections;
    for (int l = 0; l < conns.size() && !deletion_ok; ++l)
      if (conns[l].id == i) {
        swap(conns[l], conns.back());
        conns.pop_back();
        deletion_ok = true;
      }

    assert(deletion_ok);
    computeChannelThroughput(sol(ts, channel));

    channel = sol(ts, channel).parent;
  } while (channel != -1);
}

#if defined(USE_VNS_PURE)
Solution local_search(const Solution &sol20) {
  Solution best_multiple = multipleRepresentation(sol20);
  double start_of = optimal_partitioning_global(best_multiple, false);
  double best_found_before;
  do {
    best_found_before = best_multiple.throughput_;

    for (int i = 0; i < N; ++i) {
      dii best_of_i = {best_multiple.throughput_, {-2, -2}};

      Solution baseline = best_multiple;
      baseline.remove_connection(i);
      double candidate = optimal_partitioning_global(baseline);

      if (candidate > best_of_i.first)
        best_of_i = {candidate, {-1, -1}};

      for (int ts = 0; ts < sol20.slots.size(); ++ts)
        for (int c = 0; c < sol20.slots[ts].channels.size(); ++c) {
          int highest = insert_connection_up_to_root(baseline, i, ts, c);

          if (highest == -1)
            highest = c;

          candidate = optimal_partitioning_global(baseline, false);
          if (candidate > best_of_i.first)
            best_of_i = {candidate, {ts, c}};

          erase_connection_up_to_root(baseline, i, ts, c);
        }

      if (definitelyGreaterThan(best_of_i.first, best_multiple.throughput_)) {
        int ans = best_multiple.remove_connection(i);

        auto &[slot, channel] = best_of_i.second;
        if (channel < 0)
          continue;

        insert_connection_up_to_root(best_multiple, i, slot, channel);
      }
    }

    // Repeate while finding improvements
  } while (best_found_before != best_multiple.throughput_);

  double aux = optimal_partitioning_global(best_multiple);

  Solution ret{vector<TimeSlot>{TimeSlot()}};
  ret.throughput_ = best_multiple.throughput_;
  for (int t = 0; t < T; ++t)
    for (int c = 0; c < 25; ++c)
      ret(t).channels.push_back(best_multiple(t, c));

  return ret;
}

#elif defined(USE_VNS_MATH_SOLVER)

Solution local_search(const Solution &sol20) {
  Solution best_found = sol20;
  bool improved;
  do {
    improved = false;

    for (int i = 0; i < N; ++i) {
      misi mapping;

      set<int> free_connections;
      for (int j = 0; j < N; ++j)
        free_connections.insert(j);

      for (const TimeSlot &ts : sol20.slots)
        for (const Channel &ch : ts.channels)
          for (const Connection &conn : ch.connections)
            if (conn.id != i) {
              for (int o_c = 0; o_c < 45; ++o_c)
                if (overlap[ch.ix][o_c])
                  mapping[conn.id].insert(o_c);

              free_connections.erase(conn.id);
            }

      for (const int j : free_connections)
        for (int c = 0; c < 45; ++c)
          mapping[j].insert(c);

      Solution opt = vrbsp(mapping);
      if (opt > best_found) {
        // need to revert opt back to 20MHz
        // best_found = revert(opt);
        improved = true;
      }
    }

    // Repeate while finding improvements
  } while (improved);

  return best_found;
}
#endif
#endif

#if defined(USE_VNS_MATH_SOLVER) || defined(USE_VRBSP_IP) ||                   \
    defined(USE_MDVRBSP_IP)

vector<unordered_set<int>> C_b = {
    {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
        13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
    },
    {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36},
    {37, 38, 39, 40, 41, 42},
    {43, 44},
};

void couple(GRBModel *model, seila3 &y, seila2 &x) {
  for (const auto &[i, b_vars] : y)
    for (const auto &[b, mcs_vars] : y[i]) {
      vector<GRBLinExpr> y_expr_per_ts(T), x_expr_per_ts(T);
      for (const auto &[mcs, time_slots] : y[i][b])
        for (const auto &[t, _] : y[i][b][mcs])
          y_expr_per_ts[t] += y[i][b][mcs][t];

      for (const auto &[ch, time_slots] : x[i])
        for (const auto &[t, _] : x[i][ch])
          if (C_b[b].contains(ch))
            x_expr_per_ts[t] += x[i][ch][t];

      for (int t = 0; t < T; ++t) {
        string name =
            "couple" + to_string(i) + "," + to_string(b) + "," + to_string(t);
        model->addConstr(y_expr_per_ts[t] <= x_expr_per_ts[t], name);
      }
    }
}

void unique(GRBModel *model, seila2 &x) {
  for (const auto &[i, channels] : x) {
    GRBLinExpr expr = 0;
    for (const auto &[c, time_slots] : x[i])
      for (const auto &[t, var] : x[i][c])
        expr += var;

    string name = "unique" + to_string(i);
    model->addConstr(expr <= 1, name);
  }
}

void ch_overlap(GRBModel *model, seila2 &z, seila2 &x) {
  for (const auto &[i, channels] : x)
    for (const auto &[c1, _] : x[i])
      for (const auto &[t, _] : x[i][c1]) {

        GRBLinExpr expr = 0;
        for (const auto &[c2, _] : x[i])
          for (const auto &[t2, _] : x[i][c2])
            if (t == t2 && overlap[c1][c2])
              expr += x[i][c2][t2];

        string name =
            "over" + to_string(i) + "," + to_string(c1) + "," + to_string(t);
        model->addConstr(expr == z[i][c1][t], name);
      }
}

void interch(GRBModel *model, seila2 &z, seila2 &Iij) {
  for (const auto &[i, channels] : z)
    for (const auto &[c, time_slots] : z[i]) {
      for (const auto &[t, _] : z[i][c]) {

        GRBLinExpr expr = 0;
        for (const auto &[u, channels] : z) {
          if (i == u)
            continue;

          for (const auto &[c2, time_slots2] : z[u])
            for (const auto &[t2, _] : z[u][c2])
              if (t == t2 && overlap[c][c2])
                expr += AFF[u][i] * z[u][c2][t2];
        }

        string name =
            "interch" + to_string(i) + "," + to_string(c) + "," + to_string(t);
        model->addConstr(expr == Iij[i][c][t], name);
      }
    }
}

void bigG(GRBModel *model, mivar &I, seila2 &Iij, seila2 &x) {
  for (const auto &[i, time_slots] : I) {
    for (const auto &[c, time_slots] : Iij[i])
      for (const auto &[t, _] : Iij[i][c]) {
        string name =
            "bigG" + to_string(i) + "," + to_string(c) + "," + to_string(t);
        model->addConstr(I[i] >= Iij[i][c][t] - BM[i] * (1 - x[i][c][t]), name);
      }
  }
}

void bigL(GRBModel *model, mivar &I, seila2 &Iij, seila2 &x) {
  for (const auto &[i, time_slots] : I) {
    for (const auto &[c, time_slots] : Iij[i])
      for (const auto &[t, _] : Iij[i][c]) {
        string name =
            "bigL" + to_string(i) + "," + to_string(c) + "," + to_string(t);
        model->addConstr(I[i] <= Iij[i][c][t] + BM[i] * (1 - x[i][c][t]), name);
      }
  }
}

void interference(GRBModel *model, mivar &I, seila2 &x) {
  try {
    for (const auto &[i, _] : I) {
      GRBQuadExpr expr = 0;

      for (const auto &[j, _] : I) {
        if (i == j)
          continue;

        for (const auto &[c, _] : x[i]) {
          for (const auto &[t1, _] : x[i][c]) {

            for (const auto &[c2, _] : x[j]) {
              if (!overlap[c][c2])
                continue;

              for (const auto &[t2, _] : x[j][c2]) {
                if (t1 != t2)
                  continue;

                expr += AFF[j][i] * (x[i][c][t1] * x[j][c2][t1]);
              }
            }
          }
        }
      }

      string name = "interference" + to_string(i);
      model->addQConstr(I[i] == expr, name);
    }
  } catch (GRBException e) {
    cout << e.getErrorCode() << endl;
    exit(10);
  }
}

void sinr(GRBModel *model, mivar &I, seila3 &y) {
  for (const auto &[i, _] : I) {
    GRBLinExpr expr = 0;

    for (const auto &[b, mcss] : y[i]) {
      for (const auto &[m, _] : y[i][b])
        for (const auto &[t, _] : y[i][b][m])
          expr += (AFF[i][i] / SINR[m][b] - NOI) * y[i][b][m][t];
    }

    string name = "sinr" + to_string(i);
    model->addConstr(I[i] <= expr, name);
  }
}

void seila_constr(GRBModel *model, seila2 &x) {
  for (const auto &[i, channels] : x) {
    vector<vector<pair<int, GRBVar>>> vars_per_ts(T,
                                                  vector<pair<int, GRBVar>>());

    for (const auto &[c, ts] : x[i])
      for (const auto &[t, _] : x[i][c])
        vars_per_ts[t].push_back({c, x[i][c][t]});

    for (int t = 0; t < T; ++t) {
      sort(vars_per_ts[t].rbegin(), vars_per_ts[t].rend(),
           [](pair<int, GRBVar> a, pair<int, GRBVar> b) {
             return a.first < b.first;
           });

      for (int u = 0; u < vars_per_ts[t].size(); ++u)
        for (int k = u + 1; k < vars_per_ts[t].size(); ++k)
          model->addConstr(vars_per_ts[t][k].second <=
                           1 - vars_per_ts[t][u].second);
    }
  }
}

Solution vrbsp(misi &connections_channels) {
  GRBEnv env;
  GRBModel *model = new GRBModel(env);
  GRBVar *t;
  mivar I;
  seila2 x, z, Iij;
  seila3 y;

  vector<int> links;
  if (!connections_channels.empty())
    for (const auto &[l, _] : connections_channels)
      links.push_back(l);
  else
    for (int i = 0; i < N; ++i) {
      links.push_back(i);

      for (int c = 0; c < 45; ++c)
        connections_channels[i].insert(c);
    }

  // variables
  // printf("variables...\n");
  var_x(model, x, connections_channels);
  var_y(model, y, connections_channels);
  var_z(model, z, connections_channels);
  var_Iij(model, Iij, connections_channels);
  var_I(model, I, connections_channels);
  // constraints
  model->update();
  // printf("constraints...\n");
  unique(model, x);
  couple(model, y, x);
  // interference(model, I, x);
  ch_overlap(model, z, x);
  interch(model, z, Iij);
  bigG(model, I, Iij, x);
  bigL(model, I, Iij, x);
  sinr(model, I, y);
  // seila_constr(model, x);
  //
  model->update();
  model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
  model->write("model.lp");
  model->optimize();

  if (model->get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
    assert(model->get(GRB_IntAttr_SolCount) != 0);
    model->write("sol.sol");

    vector<map<int, vector<int>>> grb_sol(T);

    for (const auto &[j, channels] : x)
      for (const auto &[c, time_slots] : x[j])
        for (const auto &[t, _] : x[j][c])
          if (approximatelyEqual(x[j][c][t].get(GRB_DoubleAttr_X), 1.0))
            grb_sol[t][c].push_back(j);

    Solution ret;
    for (const map<int, vector<int>> time_slot : grb_sol) {
      TimeSlot ts;
      for (const auto &[ch, connections] : time_slot) {
        Channel channel(ChannelIdToBandwidth(ch), {});

        for (const int j : connections)
          channel.connections.push_back(Connection{j});

        ts.channels.push_back(channel);
      }

      ret.slots.push_back(ts);
    }

    computeThroughput(ret);
    return ret;
  }

  assert(false);
  return Solution();
};

// Build INTEGER PROGRAM
void var_y(GRBModel *model, seila3 &y, misi &connection_channels) {
  int cnt = 0;
  for (int t = 0; t < T; ++t) {
    for (const auto &[i, channels] : connection_channels) {
      set<int> allowed_band;
      for (const int ch : channels)
        allowed_band.insert(cToBIdx(ch));

      for (const int b : allowed_band) {
        for (int m = 0; m < 12; ++m) {
          string name = "y" + to_string(i) + "," + to_string(b) + "," +
                        to_string(m) + "," + to_string(t);
          y[i][b][m][t] = model->addVar(0.0, 1.0, DR[m][b], GRB_BINARY, name);
          cnt++;
        }
      }
    }
  }
  printf("Created %d y_vars\n", cnt);
}

void var_x(GRBModel *model, seila2 &x, misi &connection_channels) {
  int cnt = 0;
  for (int t = 0; t < T; ++t)
    for (const auto &[i, channels] : connection_channels)
      for (const int c : channels) {
        string name =
            "x" + to_string(i) + "," + to_string(c) + "," + to_string(t);
        x[i][c][t] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        cnt += 1;
      }
  printf("Created %d x_vars\n", cnt);
}

void var_t(GRBModel *model, GRBVar *t, vector<int> &links,
           misi &mapping_i_channels) {
  t = new GRBVar[T];
  for (int i = 0; i < T; ++i) {
    string name = "t" + to_string(i);
    t[i] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, name);
  }
}

void var_I(GRBModel *model, mivar &I, misi &mapping_i_channels) {
  int cnt = 0;

  for (const auto &[i, _] : mapping_i_channels) {
    string name = "I" + to_string(i);
    I[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
    cnt++;
  }

  printf("Created %d I_vars\n", cnt);
}

void var_Iij(GRBModel *model, seila2 &Iij, misi &connection_channels) {
  int cnt = 0;
  for (int t = 0; t < T; ++t) {
    for (const auto &[i, channels] : connection_channels) {
      for (const auto &c : channels) {
        string name =
            "Iij" + to_string(i) + "," + to_string(c) + "," + to_string(t);
        Iij[i][c][t] =
            model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
        cnt++;
      }
    }
  }
  printf("Created %d Iij_vars\n", cnt);
}

void var_z(GRBModel *model, seila2 &z, misi &connection_channels) {
  int cnt = 0;
  for (int t = 0; t < T; ++t)
    for (const auto &[i, channels] : connection_channels)
      for (const int c : channels) {
        string name =
            "z" + to_string(i) + "," + to_string(c) + "," + to_string(t);
        z[i][c][t] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        cnt++;
      }

  printf("Created %d z_vars\n", cnt);
}

#endif

int bToIdx(int b) {
  if (b == 20)
    return 0;
  else if (b == 40)
    return 1;
  else if (b == 80)
    return 2;

  return 3;
}

bool GreaterOrEqual(double a, double b) {
  return definitelyGreaterThan(a, b) || approximatelyEqual(a, b);
}

int necessaryBw(double dr) {
  for (int i = 0; i < 4; ++i) {
    if (definitelyGreaterThan(DR[11][i], dr))
      return i;
  }

  return -1;
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

int idxToB(int idx) {
  if (idx == 0)
    return 20;
  else if (idx == 1)
    return 40;
  else if (idx == 2)
    return 80;

  return 160;
}

int ChannelIdToBandwidth(int c) {
  if (c <= 24)
    return 20;
  else if (c <= 36)
    return 40;
  else if (c <= 42)
    return 80;
  else
    return 160;
}

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

bool canTransmitUsingChannel(int i, int c) {
  return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]) ||
         approximatelyEqual(DR[11][cToBIdx(c)], GMM[i]);
}

bool canTransmitUsingBandwidth(int i, int b, int m) {
  return definitelyGreaterThan(DR[m][b], GMM[i]) ||
         approximatelyEqual(DR[m][b], GMM[i]);
}

vector<string> gurobi_sol(const Solution &sol) {
  vector<string> ans;
  for (int ts_ix = 0; const TimeSlot &ts : sol.slots) {
    for (const Channel &ch : ts.channels)
      for (const Connection &conn : ch.connections)
        ans.push_back("x" + to_string(conn.id) + "," + to_string(ch.ix) + "," +
                      to_string(ts_ix));

    ts_ix++;
  }

  return ans;
}
