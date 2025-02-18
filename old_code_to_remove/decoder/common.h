#ifndef COMMON
#define COMMON

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

const int C = 45;
const int MAX_SPECTRUM = 5;
const int MAX_CHANNELS = 45;
const int MAX_SLOTS = 2048;

extern int N, T, n_spectrums;
extern double alfa, NOI, powerSender;
extern double receivers[MAX_CONN][2], senders[MAX_CONN][2];
extern double AFF[MAX_CONN][MAX_CONN];
extern double DM[MAX_CONN][MAX_CONN]; //, interferenceMatrix[MAX_CONN][MAX_CONN];
extern vector<vector<double>> DR, SINR, B;
extern vector<int> spectrum_size;

extern double maximumTime;

extern vector<double> GMM, BM;

extern random_device rd;
extern default_random_engine whatever;

extern int overlap[45][45];

double distance(double, double, double, double);

bool approximatelyEqual(double, double, double epsilon = EPS);

bool essentiallyEqual(double, double, double epsilon = EPS);

bool definitelyGreaterThan(double, double, double epsilon = EPS);

bool definitelyLessThan(double, double, double epsilon = EPS);

bool GreaterOrEqual(double, double);

void convertTableToMW(const vector<vector<double>> &, vector<vector<double>> &);

void distanceAndInterference();

double convertDBMToMW(double _value);

double gammaToBeta(double gmm, int bw);

int cToBIdx(int c);

int bToIdx(int b);

int idxToB(int idx);

int necessaryBw(double dr);

int avgNts();

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

    Channel(double throughput, double interference, double violation, int bandwidth,
            const vector<Connection> connections);
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

    double decode(std::vector<double> &) const;

    TimeSlot &operator()(int t);
    Channel &operator()(int t, int c);
    Connection &operator()(int t, int c, int l);
};

void read_data();

#endif
