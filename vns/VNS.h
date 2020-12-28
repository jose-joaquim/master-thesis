#ifndef VNS_H
#define VNS_H

#include "MTRand.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <numeric>
#include <random>
#include <thread>
#include <vector>

const int MAX_CONN = 2048;
const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-9;

using ii = std::pair<int, int>;
using dd = std::pair<double, double>;

static MTRand rng;
static int n_connections, n_spectrums;
static double alfa, noise, powerSender;
static std::vector<std::vector<double>> dataRates[4], SINR[4];
static double receivers[MAX_CONN][2], senders[MAX_CONN][2];
static double affectance[MAX_CONN][MAX_CONN];
static double distanceMatrix[MAX_CONN][MAX_CONN], interferenceMatrix[MAX_CONN][MAX_CONN];

static const int MAX_SPECTRUM = 4;
static const int MAX_CHANNELS = 45;
static double startTime;
static clock_t maximumTime;

#ifdef MDVRBSP
static std::vector<double> gamma;
#endif

class Solution {};

void vns(void);

void constructive_heuristic(void);

void local_search(void);

void init(void);

bool double_equals(double a, double b, double epsilon) { return std::abs(a - b) < epsilon; }

bool stop(void) { return (((double)(clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime; }

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
    return hypot((X_si - X_ri), (Y_si - Y_ri));
}

void convertTableToMW(std::vector<std::vector<double>> (&_SINR)[4],
                      std::vector<std::vector<double>> (&__SINR)[4]) {}

void distanceAndInterference() {
    for (int i = 0; i < n_connections; i++) {
        double X_si = receivers[i][X_c];
        double Y_si = receivers[i][Y_c];

        for (int j = 0; j < n_connections; j++) {

            double X_rj = senders[j][X_c];
            double Y_rj = senders[j][Y_c];

            double dist = distance(X_si, Y_si, X_rj, Y_rj);

            distanceMatrix[i][j] = dist;

            if (i == j) {
                interferenceMatrix[i][j] = 0.0;
            } else {
                double value = (dist != 0.0) ? powerSender / pow(dist, alfa) : 1e9;
                interferenceMatrix[i][j] = value;
            }
        }
    }
}

void read_data(const int type) {
    scanf("%d %lf %lf %lf %d", &n_connections, &alfa, &noise, &powerSender, &n_spectrums);

    for (int i = 0; i < n_spectrums; i++) {
        int s;
        scanf("%d", &s);
        // TODO: write the spectrum size...
    }

    if (noise != 0) {
        // TODO: convert DBM to MW
    }

    // Read receivers
    for (int i = 0; i < n_connections; i++) {
        double a, b;
        scanf("%lf %lf", &a, &b);
        receivers[i][0] = a;
        receivers[i][1] = b;
    }

    // Read Senders
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
        gamma.emplace_back(bt);
#endif
    }

    for (int i = 0; i < 12; i++) {
        std::vector<double> aux;
        for (int j = 0; j < 4; j++) {
            double a;
            scanf("%lf", &a);
            aux.emplace_back(a);
        }
        dataRates[i].emplace_back(aux);
    }

    for (int i = 0; i < 12; i++) {
        std::vector<double> aux;
        for (int j = 0; j < 4; j++) {
            double a;
            scanf("%lf", &a);
            aux.emplace_back(a);
        }
        SINR[i].emplace_back(aux);
    }

    convertTableToMW(SINR, SINR);
    distanceAndInterference();

    for (int i = 0; i < n_connections; i++) {
        for (int j = 0; j < n_connections; j++) {
            // TODO: affectance
        }
    }

    for (int i = 0; i < n_connections; i++) {
        for (int j = 0; j < 4; j++) {
            // TODO: beta
        }
    }
}

#endif
