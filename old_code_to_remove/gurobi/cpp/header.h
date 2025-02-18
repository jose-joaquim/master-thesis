#ifndef VNS_H
#define VNS_H

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
const double EPS = 1e-6;

using namespace std;
using namespace std::chrono;

using ii = pair<int, int>;
using dd = pair<double, double>;
using ti3 = tuple<int, int, int>;

static const int C = 45;
static int N, T, n_spectrums;
static double alfa, NOI, powerSender;
static double receivers[MAX_CONN][2], senders[MAX_CONN][2];
static double AFF[MAX_CONN][MAX_CONN];
static double DM[MAX_CONN]
                [MAX_CONN];  //, interferenceMatrix[MAX_CONN][MAX_CONN];
static vector<vector<double>> DR, SINR, B;
static vector<int> spectrum_size;

static const int MAX_SPECTRUM = 5;
static const int MAX_CHANNELS = 45;
static const int MAX_SLOTS = 2048;
static double maximumTime;

static vector<double> GMM, BM;

random_device rd;
auto whatever = default_random_engine{rd()};

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

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
  return hypot((X_si - X_ri), (Y_si - Y_ri));
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

bool GreaterOrEqual(double a, double b) {
  return definitelyGreaterThan(a, b) || approximatelyEqual(a, b);
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
  int mcs = 0;

  while (mcs < 12 && definitelyGreaterThan(gmm, DR[mcs][bw])) mcs++;

  if (mcs < 12) return SINR[mcs][bw];

  return SINR[11][3] + 1.0;
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

int idxToB(int idx) {
  if (idx == 0)
    return 20;
  else if (idx == 1)
    return 40;
  else if (idx == 2)
    return 80;

  return 160;
}

int necessaryBw(double dr) {
  for (int i = 0; i < 4; ++i) {
    if (definitelyGreaterThan(DR[11][i], dr)) return i;
  }

  return -1;
}

int avgNts() {
  double ret = 0;
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < 4; ++j) {
      if (definitelyGreaterThan(DR[11][j], GMM[i])) {
        if (j == 0)
          ret += 20;
        else if (j == 1)
          ret += 40;
        else if (j == 2)
          ret += 80;
        else
          ret += 160;
      }
    }
  }
  return ceil(ret / 500.0);
}

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

struct TimeSlot {
  vector<Channel> channels;
  double interference;
  double throughput;

  TimeSlot(const vector<Channel> &);
  TimeSlot();

  Channel &operator()(int c);
  Connection &operator()(int c, int l);
};

class Solution {
 public:
  vector<TimeSlot> slots;
  double throughput;
  double violation;

  Solution(const vector<Channel> sp, double tot, bool flag);
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
  Channel &operator()(int t, int c);
  Connection &operator()(int t, int c, int l);
};

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

TimeSlot::TimeSlot(const vector<Channel> &sp) : TimeSlot() { channels = sp; }

TimeSlot::TimeSlot() {
  interference = 0.0;
  throughput = 0.0;
  channels = vector<Channel>();
}

Channel &TimeSlot::operator()(int c) { return channels[c]; }

Connection &TimeSlot::operator()(int c, int l) {
  return channels[c].connections[l];
}

Solution::Solution(const vector<Channel> sp, double tot, bool flag)
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

Channel &Solution::operator()(int t, int c) { return slots[t].channels[c]; }

Connection &Solution::operator()(int t, int c, int l) {
  return slots[t].channels[c].connections[l];
}

bool SINRfeas(const Connection &con, int bw) {
  return GreaterOrEqual(con.SINR, B[con.id][bToIdx(bw)]);
}

double computeConnectionThroughput(Connection &i, int bw) {
  if (bw == 0) return 0.0;

  if (approximatelyEqual(i.interference, 0.0, EPS)) {
    i.SINR = SINR[11][bToIdx(bw)];
    i.throughput = DR[11][bToIdx(bw)];
  } else {
    i.SINR = AFF[i.id][i.id] / (i.interference + NOI);

    int mcs = 11;
    while (mcs >= 0 && definitelyLessThan(i.SINR, SINR[mcs][bToIdx(bw)])) mcs--;

    if (mcs < 0) {
      return -1;
    } else
      i.throughput = DR[mcs][bToIdx(bw)];
  }

  return i.throughput;
}

#ifdef MDVRBSP
optional<Channel> insertInChannel(Channel new_c, int i_idx) {
  Connection i(i_idx, 0.0, 0.0, DM[i_idx][i_idx]);

  for (auto &j : new_c.connections) {
    j.interference += AFF[i_idx][j.id];
    i.interference += AFF[j.id][i_idx];
  }

  new_c.connections.emplace_back(i);
  for (auto &i : new_c.connections) {
    double ans = computeConnectionThroughput(i, new_c.bandwidth);
    bool feas = SINRfeas(i, new_c.bandwidth);
    // if (feas) {
    //   if (definitelyGreaterThan(GMM[i.id], ans)) {
    //     printf("%d %.4lf %.4lf %.4lf %.4lf %d\n", feas, i.SINR,
    //            B[i.id][bToIdx(new_c.bandwidth)], i.throughput, GMM[i.id],
    //            new_c.bandwidth);
    //     assert(false);
    //   }
    // }
    if (!SINRfeas(i, new_c.bandwidth)) return nullopt;

    new_c.throughput += i.throughput;
  }

  return new_c;
}

optional<vector<Channel>> split(Channel toSplit, vector<double> *GMA) {
  vector<Channel> ret(2, Channel(toSplit.bandwidth / 2));

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

    if (!bestChannel) return nullopt;

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
  if (!ch) return nullopt;

  // for (const Connection &conn : ch->connections) {
  //   cout << conn.id << " " << conn.SINR << " " << conn.interference << " "
  //        << B[conn.id][bToIdx(ch->bandwidth)] << endl;
  //
  //   assert(definitelyGreaterThan(conn.throughput, GMM[conn.id]) ||
  //          approximatelyEqual(conn.throughput, GMM[conn.id]));
  //   assert(
  //       definitelyGreaterThan(conn.SINR, B[conn.id][bToIdx(ch->bandwidth)])
  //       || approximatelyEqual(conn.SINR,
  //       B[conn.id][bToIdx(ch->bandwidth)]));
  // }

  for (Connection &conn : ch->connections) {
    if (definitelyLessThan(conn.SINR, B[conn.id][bToIdx(ch->bandwidth)]))
      return nullopt;
    // assert(definitelyGreaterThan(conn.SINR,
    // B[conn.id][bToIdx(ch->bandwidth)]));
  }

  return ch;
}

bool can_split(const Channel &ch) {
  if (ch.bandwidth <= 20) return false;

  for (const auto &conn : ch.connections)
    if (definitelyLessThan(DR[11][bToIdx(ch.bandwidth / 2)], GMM[conn.id]))
      return false;

  return true;
}

Solution CH() {
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
      if (conn < t) continue;
      for (int c = 0; c < ret(t).channels.size(); c++) {
        auto rst = try_insert(conn, ret(t, c));

        if (rst) {
          ret(t, c) = *rst;
          // for (int t = 0; const TimeSlot &ts : ret.slots) {
          //   for (const Channel &ch : ts.channels) {
          //     for (const Connection &conn : ch.connections) {
          //       if (!SINRfeas(conn, ch.bandwidth)) {
          //         printf("%.4lf %.4lf %.4lf %.4lf\n", conn.SINR,
          //                B[conn.id][bToIdx(ch.bandwidth)], conn.throughput,
          //                GMM[conn.id]);
          //       }
          //
          //       assert(SINRfeas(conn, ch.bandwidth));
          //     }
          //   }
          // }

          goto NEXT_CONN;
        }

        if (can_split(ret(t, c))) {
          auto ch_split = split(ret(t, c), &GMM);

          if (!ch_split) continue;

          for (int i = 0; i < ch_split->size(); ++i) {
            auto cp_channel = try_insert(conn, (*ch_split)[i]);

            if (cp_channel) {
              swap(ret(t, c), ret(t).channels.back());
              ret(t).channels.pop_back();
              ret(t).channels.emplace_back(*cp_channel);
              ret(t).channels.emplace_back((*ch_split)[i == 0]);
              // for (int t = 0; const TimeSlot &ts : ret.slots) {
              //   for (const Channel &ch : ts.channels) {
              //     for (const Connection &conn : ch.connections) {
              //       if (!SINRfeas(conn, ch.bandwidth)) {
              //         printf("%.4lf %.4lf %.4lf %.4lf\n", conn.SINR,
              //                B[conn.id][bToIdx(ch.bandwidth)],
              //                conn.throughput, GMM[conn.id]);
              //       }
              //       assert(SINRfeas(conn, ch.bandwidth));
              //     }
              //   }
              // }
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

  for (int t = 0; const TimeSlot &ts : ret.slots) {
    for (const Channel &ch : ts.channels) {
      for (const Connection &conn : ch.connections) {
        assert(SINRfeas(conn, ch.bandwidth));
      }
    }
  }

  return ret;
}

vector<tuple<int, int, int>> gurobi_sol(const Solution &sol) {
  vector<vector<int>> C_b = {
      {
          0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
          13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
      },
      {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36},
      {37, 38, 39, 40, 41, 42},
      {43, 44},
  };

  vector<tuple<int, int, int>> ans;
  set<int> allow_tmp;
  for (int i = 0; i < 45; ++i) allow_tmp.insert(i);

  vector<set<int>> channels_free(sol.slots.size(), allow_tmp);

  for (int ts_idx = 0; const TimeSlot &ts : sol.slots) {
    for (const Channel &ch : ts.channels) {
      if (ch.connections.empty()) continue;

      bool inserted = false;
      for (const int candidate : C_b[bToIdx(ch.bandwidth)]) {
        if (channels_free[ts_idx].contains(candidate)) {
          // printf("CHANNEL (%d, %d, %d):\n", candidate, i, ch.bandwidth);
          for (const Connection &conn : ch.connections) {
            // cout << "     " << conn.id << " " << conn.interference <<
            // endl;
            ans.push_back({conn.id, candidate, ts_idx});
            inserted = true;
          }
          // printf("end\n");

          for (int j = 0; j < 45; ++j)
            if (overlap[candidate][j]) channels_free[ts_idx].erase(j);

          break;
        }
      }
      assert(inserted);
    }
    ts_idx++;
  }

  return ans;
}
#endif

void read_data() {
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

  for (int i = 0; i < N; i++) {
    double bt;
    scanf("%lf", &bt);
#ifdef MDVRBSP
    GMM.emplace_back(bt);
#else
    GMM.emplace_back(0.0);
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

  BM.assign(N, 0.0);
  B.assign(N, vector<double>(4, 0));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; ++j) {
      if (i != j) BM[i] += AFF[j][i];

      if (i != j) {
        int bwj = necessaryBw(GMM[j]);
        int bwi = necessaryBw(GMM[i]);

        if (bwj >= bwi) BM[i] += AFF[j][i] + 0.0005;
      }
    }

    for (int j = 0; j < 4; j++) {
      double value = gammaToBeta(GMM[i], j);
      assert(definitelyGreaterThan(value, 0.0));
      B[i][j] = value;
    }
  }
}

#endif
