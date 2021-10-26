//
// Created by José Joaquim on 03/03/20.
//

#ifndef BRKGA_FF_BEST_HEURISTICDECODER_H
#define BRKGA_FF_BEST_HEURISTICDECODER_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include <chrono>
#include <thread>

#include "MTRand.h"

using namespace std;

const int MAX_CONN = 2048;

extern int nConnections, nSpectrums;
extern double senders[MAX_CONN][2], receivers[MAX_CONN][2];
extern pair<int, int> zeroChannel;
extern double powerSender, alfa, noise, ttm;
extern MTRand rng;
extern const int MAX_SPECTRUM, MAX_CHANNELS;

using ii = pair<int, int>;

struct Connection {
    int id;
    double throughput;
    double interference;
    double SINR;
    double distanceSR;

    Connection(int id, double throughput, double interference, double distanceSR);

    Connection(int id);

    bool operator<(const Connection &other) const { return distanceSR < other.distanceSR; }

    bool operator>(const Connection &other) const { return !operator<(other); }
};

struct Channel {
    double throughput; // TODO: So far, I do not have any use for this variable.
    double interference;
    int bandwidth;
    vector<Connection> connections;

    bool operator<(const Channel &other) const { return bandwidth < other.bandwidth; }

    bool operator>(const Channel &other) const { return !operator<(other); }

    Channel(double throughput, double interference, int bandwidth,
            const vector<Connection> &connections);

    Channel(int bandwidth);
};

struct Spectrum {
    int maxFrequency;
    int usedFrequency;
    vector<Channel> channels;

    Spectrum(int maxFrequency, int usedFrequency, const vector<Channel> &channels);
};

class Solution {
  public:
    vector<Spectrum> spectrums;
    double totalThroughput;

    bool throughputFlag; // FIXME: conflict with the constructors

    Solution(const vector<Spectrum> &spectrums, double totalThroughput, bool throughputFlag);

    Solution();

    void printSolution(FILE *file = nullptr);

    // AVAILABLE ONLY IN THE META-HEURISTICS FILES
    double decode(vector<double> variables) const;

    friend bool operator>(const Solution &o1, const Solution &o2);

    friend bool operator<(const Solution &o1, const Solution &o2);
};

bool checkOne(const Solution &s);

bool checkTwo(const Solution &s);

int bwIdx(int bw);

void loadData();

void rawInsert(Solution &sol, int conn, ii where);

void rawRemove(Solution &sol, int conn, ii where);

Channel insertInChannel(Channel newChannel, int conn);

Channel deleteFromChannel(const Channel &channel, int conn);

double computeThroughput(Solution &curr, bool force = false);

bool double_equals(double a, double b, double epsilon = 0.000000001);

void computeChannelsThroughput(vector<Channel> &channels);

double computeConnectionThroughput(Connection &conn, int bandWidth, bool force = false);

Solution createSolution();

int computeConnectionMCS(const Connection &conn, int bandwidth);

#endif // BRKGA_FF_BEST_HEURISTICDECODER_H
