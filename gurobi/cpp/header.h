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
const double EPS = 1e-9;

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

static high_resolution_clock::time_point start;

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
    distanceSR = DM[id][id];
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

  bool operator<(const Solution &o1) const { return violation < o1.violation; }

  bool operator>(const Solution &o1) const { return !operator<(o1); }

  Channel &operator()(int t, int s, int c) {
    return slots[t].spectrums[s].channels[c];
  }
};

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

double computeConnectionThroughput(Connection &conn, int bandwidth) {
  if (bandwidth == 0) return 0.0;

  int mcs = -1;
  int maxDataRate = 11;

  if (approximatelyEqual(conn.interference, 0.0, EPS)) {
    mcs = maxDataRate;
    conn.throughput = DR[mcs][bwIdx(bandwidth)];
  } else {
    conn.SINR = AFF[conn.id][conn.id] / (conn.interference + NOI);

    mcs = 11;
    while (mcs >= 0 && conn.SINR < SINR[mcs][bwIdx(bandwidth)]) mcs--;

    if (mcs < 0)
      conn.throughput = 0.0;
    else
      conn.throughput = DR[mcs][bwIdx(bandwidth)];
  }

  return conn.throughput;
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

double computeThroughput(Solution &curr, bool force = false) {
  double OF = 0.0;

  for (int t = 0; t < curr.slots.size(); t++) {
    for (int s = 0; s < curr.slots[t].spectrums.size(); s++) {
      for (int c = 0; c < curr.slots[t].spectrums[s].channels.size(); c++) {
        auto &chan = curr.slots[t].spectrums[s].channels[c];
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

Channel insertLink(Channel newChannel, int idConn) {  // TODO: remover macro
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
      Channel inserted = insertLink(ret[c], connToInsert[i].id);
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

bool is_feasible(const Solution &ret, bool retOption = true) {
  for (const TimeSlot &ts : ret.slots) {
    for (const Spectrum &sp : ts.spectrums) {
      for (const Channel &ch : sp.channels) {
        for (const Connection &conn : ch.connections) {
          if (definitelyLessThan(conn.throughput, GMM[conn.id])) {
            printf("ops... conn %d %lf %lf\n", conn.id, conn.throughput,
                   GMM[conn.id]);
            return false;
          }
        }
      }
    }
  }

  return true;
}

// int feas_recur(const Channel &c) {
//     printf("recur bw %d |c| %d\n", c.bandwidth, (int)c.connections.size());
//     int ret = c.bandwidth;
//     if (!is_feasible(c)) {
//         int bw1 = c.bandwidth, bw2 = c.bandwidth > 20 ? c.bandwidth / 2 :
//         c.bandwidth; Channel c1(bw1), c2(bw2);
//
//         vector<int> links1, links2;
//         for (int i = 0; i < c.connections.size(); ++i) {
//             int a = necessaryBw(GMM[c(i)]);
//             int b = idxToB(a);
//             cout << a << " " << c(i) << " " << b << endl;
//             if (b < c.bandwidth)
//                 links2.emplace_back(c(i));
//             else links1.emplace_back(c(i));
//         }
//
//         printf("links1 eh %d links2 eh %d\n", (int)links1.size(),
//         (int)links2.size()); for (int i = 0; i < links1.size(); ++i)
//             c1 = insertLink(c1, links1[i]);
//
//         for (int i = 0; i < links2.size(); ++i)
//             c2 = insertLink(c2, links2[i]);
//
//         ret = feas_recur(c1);
//         ret += feas_recur(c2);
//     }
//
//     printf("ret eh %d\n", ret);
//     return ret;
// }
//
// int aux() {
//     int ans = 0;
//     map<int, vector<int>> lPerBw;
//
//     int maxBw = -1;
//     for (int i = 0; i < N; ++i)
//         maxBw = max(maxBw, idxToB(necessaryBw(GMM[i])));
//
//     Channel c(maxBw);
//     for (int i = 0; i < N; ++i) {
//         c = insertLink(c, i);
//         // int bw = idxToB(necessaryBw(GMM[i]));
//         // lPerBw[bw].emplace_back(i);
//     }
//
//     ans += feas_recur(c);
//
//     // for (auto &e : lPerBw) {
//     //     printf("channel of %d: ", e.first);
//     //     Channel c(e.first);
//     //     for (int i = 0; i < e.second.size(); ++i) {
//     //         c = insertLink(c, e.second[i]);
//     //     }
//     //
//     //     ans += feas_recur(c);
//     // }
//
//     return ceil(ans * 1.0 / 500.0);
//
// }

bool try_insert(int conn_id, Channel &ch) {
  ch = insertLink(ch, conn_id);

  for (Connection &conn : ch.connections)
    if (definitelyLessThan(conn.throughput, GMM[conn.id])) {
      return false;
    }

  return true;
}

bool can_split(const Channel &ch) {
  if (ch.bandwidth <= 20) return false;

  for (const auto &conn : ch.connections)
    if (definitelyLessThan(DR[11][bToIdx(ch.bandwidth)], GMM[conn.id]))
      return false;

  Channel aux(ch.bandwidth / 2);
  for (auto &conn : ch.connections)
    if (!try_insert(conn.id, aux)) return false;

  return true;
}

#ifdef MDVRBSP
int CH() {
  TimeSlot dummy_ts(init_conf);
  Solution ret({dummy_ts});

  vector<int> links;
  for (int i = 0; i < N; ++i) links.emplace_back(i);

  shuffle(links.begin(), links.end(), whatever);

  for (int conn : links) {
    bool success = false;
    for (int t = 0; t < ret.slots.size(); t++) {
      for (int s = 0; s < ret.slots[t].spectrums.size(); s++) {
        if (success) break;

        for (int c = 0; c < ret.slots[t].spectrums[s].channels.size(); c++) {
          if (success) break;

          Channel cp_channel = ret(t, s, c);
          success = try_insert(conn, cp_channel);

          if (success) {
            ret(t, s, c) = cp_channel;
            break;
          }

          if (can_split(ret(t, s, c))) {
            assert(approximatelyEqual(ret(t, s, c).violation, 0.0));
            vector<Channel> ch_split = split(ret(t, s, c));

            assert(approximatelyEqual(ch_split[0].violation, 0.0));
            assert(approximatelyEqual(ch_split[1].violation, 0.0));

            cp_channel = ch_split[0];
            success = try_insert(conn, cp_channel);

            if (success) {
              swap(ret.slots[t].spectrums[s].channels[c],
                   ret.slots[t].spectrums[s].channels.back());
              ret.slots[t].spectrums[s].channels.pop_back();
              ret.slots[t].spectrums[s].channels.emplace_back(cp_channel);
              ret.slots[t].spectrums[s].channels.emplace_back(ch_split[1]);
              break;
            }

            cp_channel = ch_split[1];
            success = try_insert(conn, cp_channel);

            if (success) {
              swap(ret.slots[t].spectrums[s].channels[c],
                   ret.slots[t].spectrums[s].channels.back());
              ret.slots[t].spectrums[s].channels.pop_back();
              ret.slots[t].spectrums[s].channels.emplace_back(cp_channel);
              ret.slots[t].spectrums[s].channels.emplace_back(ch_split[0]);
              break;
            }
          }
        }
      }
    }

    if (!success) {
      TimeSlot new_ts(dummy_ts);
      for (int s = 0; s < new_ts.spectrums.size() && !success; s++) {
        for (int c = 0; c < new_ts.spectrums[s].channels.size() && !success;
             c++) {
          int bw = new_ts.spectrums[s].channels[c].bandwidth;
          if (bw >= 160) {
            Connection new_conn(conn);
            new_conn.SINR = 1000000007;
            new_conn.throughput = DR[11][bwIdx(bw)];
            new_ts.spectrums[s].channels[c].connections.push_back(new_conn);
            ret.slots.emplace_back(new_ts);
            success = true;
            break;
          }
        }
      }
    }
    computeThroughput(ret);
    assert(is_feasible(ret));
  }

  computeThroughput(ret);
  return (int)ret.slots.size();
}
#endif

void read_data() {
  scanf("%d %lf %lf %lf %d", &N, &alfa, &NOI, &powerSender, &n_spectrums);

  for (int i = 0; i < n_spectrums; i++) {
    int s;
    scanf("%d", &s);
    spectrum_size.emplace_back(s);
    init_conf.emplace_back(s, 0, vector<Channel>());
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
    GMM.emplace_back(bt);
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

  initTimeSlot();
  printf("%d %d\n", avgNts(), CH());
  // T = avgNts();
#ifdef MDVRBSP
  T = CH();
  cout << T << endl;
#endif
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
      B[i][j] = value;
    }
  }
}

#endif
