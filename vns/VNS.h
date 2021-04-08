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

using ii = std::pair<int, int>;
using dd = std::pair<double, double>;
using ti3 = std::tuple<int, int, int>;

static MTRand rng;
static int n_connections, n_spectrums, n_time_slots;
static double alfa, noise, powerSender;
static double receivers[MAX_CONN][2], senders[MAX_CONN][2];
static double affectance[MAX_CONN][MAX_CONN];
static double distanceMatrix[MAX_CONN][MAX_CONN], interferenceMatrix[MAX_CONN][MAX_CONN];
static std::vector<std::vector<double>> dataRates, SINR, beta;
static std::vector<int> spectrum_size;

static const int MAX_SPECTRUM = 4;
static const int MAX_CHANNELS = 45;
static const int MAX_SLOTS = 2048; // TODO: This is a high value
static double startTime;
static clock_t maximumTime;

static int parent[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
static int child[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS][2];
static double chanThroughput[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
static bool inSolution[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];
static double sc_opt[MAX_SLOTS][MAX_SPECTRUM][MAX_CHANNELS];

#ifdef MDVRBSP
static std::vector<double> gma;
#endif

std::random_device rd;
auto whatever = std::default_random_engine{rd()};

struct Connection {
    int id;
    double throughput;
    double interference;
    double SINR;
    double distanceSR;

    Connection(int id, double throughput, double interference, double distanceSR)
        : id(id), throughput(throughput), interference(interference), distanceSR(distanceSR) {}

    Connection(int id) : id(id) {
        throughput = 0.0;
        interference = 0.0;
        SINR = 0.0;
        distanceSR = distanceMatrix[id][id];
    }

    bool operator<(const Connection &other) const { return distanceSR < other.distanceSR; }

    bool operator>(const Connection &other) const { return !operator<(other); }
};

struct Channel {
    double throughput; // TODO: So far, I do not have any use for this variable.
    double interference;
    int bandwidth;
    std::vector<Connection> connections;

    bool operator<(const Channel &other) const { return bandwidth < other.bandwidth; }

    bool operator>(const Channel &other) const { return !operator<(other); }

    Channel(double throughput, double interference, int bandwidth,
            const std::vector<Connection> &connections)
        : throughput(throughput), interference(interference), bandwidth(bandwidth),
          connections(connections) {}

    Channel(int bandwidth) : bandwidth(bandwidth) {
        interference = 0.0;
        throughput = 0.0;
        connections = std::vector<Connection>();
    }

    Channel() : Channel(0.0) {}
};

struct Spectrum {
    int maxFrequency;
    int usedFrequency;
    double interference;
    std::vector<Channel> channels;

    Spectrum(int maxFrequency, int usedFrequency, const std::vector<Channel> &channels)
        : maxFrequency(maxFrequency), usedFrequency(usedFrequency), channels(channels) {}

    Spectrum() {
        maxFrequency = 0;
        usedFrequency = 0;
        interference = 0.0;
        channels = std::vector<Channel>();
    }
};

struct TimeSlot {
    std::vector<Spectrum> spectrums;
    double interference;
    double throughput;

    TimeSlot(const std::vector<Spectrum> &sp) : spectrums(sp) {}
    TimeSlot() {
        interference = 0.0;
        throughput = 0.0;
        spectrums = std::vector<Spectrum>();
    }
};

class Solution {
  public:
    std::vector<TimeSlot> slots;
    std::vector<int> unscheduled;
    double throughput;
    double violation;
    bool throughput_flag;

    Solution(const std::vector<Spectrum> &sp, double tot, bool flag)
        : throughput(tot), throughput_flag(flag) {
        slots.push_back(TimeSlot(sp));
    }

    Solution(const vector<TimeSlot> &ts) : slots(ts) {}

    Solution() {}

#ifdef MDVRBSP
    bool operator<(const Solution &o1) const { return violation < o1.violation; }

    bool operator>(const Solution &o1) const { return !operator<(o1); }
#else
    bool operator<(const Solution &o1) const { return throughput < o1.throughput; }

    bool operator>(const Solution &o1) const { return !operator<(o1); }
#endif
};

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

#ifdef MDVRBSP
double compute_violation(Solution &sol) {
    sol.violation = 0.0;
    for (const TimeSlot &ts : sol.slots) {
        for (const Spectrum &sp : ts.spectrums) {
            for (const Channel &ch : sp.channels) {
                for (const Connection &conn : ch.connections) {
                    double res = gma[conn.id] - conn.throughput;
                    sol.violation = max(sol.violation, gma[conn.id] - conn.throughput);
                }
            }
        }
    }
    return sol.violation;
}
#endif

bool approximatelyEqual(double a, double b, double epsilon = EPS)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentiallyEqual(double a, double b, double epsilon = EPS)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyGreaterThan(double a, double b, double epsilon = EPS)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitelyLessThan(double a, double b, double epsilon = EPS)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool stop(void) { return (((double)(clock() - startTime)) / CLOCKS_PER_SEC) >= maximumTime; }

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
    return hypot((X_si - X_ri), (Y_si - Y_ri));
}

void convertTableToMW(const std::vector<std::vector<double>> &_SINR,
                      std::vector<std::vector<double>> &__SINR) {
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

double convertDBMToMW(double _value) {
    double result = 0.0, b;

    b = _value / 10.0;     // dBm dividido por 10
    result = pow(10.0, b); // Converte de DBm para mW

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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
static std::vector<Spectrum> init_conf;
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
                sp.channels.emplace_back(0.0, 0.0, bw, vector<Connection>());
            }
        }
        assert(sp.maxFrequency - sp.usedFrequency >= 0);
    }
}

void read_data() {
    scanf("%d %lf %lf %lf %d", &n_connections, &alfa, &noise, &powerSender, &n_spectrums);

    for (int i = 0; i < n_spectrums; i++) {
        int s;
        scanf("%d", &s);
        spectrum_size.emplace_back(s);
        init_conf.emplace_back(s, 0, std::vector<Channel>());
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
        gma.emplace_back(bt);
#endif
    }

    dataRates.assign(12, std::vector<double>(4, 0));
    for (int i = 0; i < 12; i++) {
        for (int j = 0; j < 4; j++) {
            double a;
            scanf("%lf", &a);
            dataRates[i][j] = a;
        }
    }

    SINR.assign(12, std::vector<double>(4, 0));
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

#ifdef MDVRBSP
    beta.assign(n_connections, std::vector<double>(4, 0));
    for (int i = 0; i < n_connections; i++) {
        for (int j = 0; j < 4; j++) {
            double value = gammaToBeta(gma[i], j);
            beta[i][j] = value;
        }
    }
#endif
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

void print_solution_to_gurobi(const Solution &sol, FILE *file = nullptr) {
    if (file == nullptr) {
        file = fopen("./solution_gurobi.txt", "w");
    }

    fprintf(file, "%d\n", int(sol.slots.size()));
    
    for (const TimeSlot &ts : sol.slots) {
        vector<string> lines;
        for (const Spectrum &sp : ts.spectrums) {
            for (const Channel &ch : sp.channels) {
                if (ch.connections.empty())
                    continue;
                
                string line = to_string(ch.bandwidth) + " ";
                
                for (int i = 0; i < ch.connections.size() - 1; i++) {
                    line += (to_string(ch.connections[i].id) + " ");
                }

                line += to_string(ch.connections.back().id);
                lines.emplace_back(line + "\n");
            }
        }

        fprintf(file, "%d\n", int(lines.size()));
        for (const string &line : lines) {
            fprintf(file, "%s", line.c_str());
        }
    }

    
    
    fclose(file);
}

void print_solution_to_file(const Solution &sol, FILE *file = nullptr) {
    if (file == nullptr) {
        file = fopen("./solution.txt", "w");
    }

    fprintf(file, "%lf %d\n", sol.throughput, int(sol.slots.size()));
    for (int t = 0; t < sol.slots.size(); t++) {
        fprintf(file, "%d\n", int(sol.slots[t].spectrums.size()));
        for (int s = 0; s < sol.slots[t].spectrums.size(); s++) {
            fprintf(file, "%d\n", int(sol.slots[t].spectrums[s].channels.size()));
            for (int c = 0; c < int(sol.slots[t].spectrums[s].channels.size()); c++) {
                fprintf(file, "%d %d", sol.slots[t].spectrums[s].channels[c].bandwidth,
                        int(sol.slots[t].spectrums[s].channels[c].connections.size()));
                for (int i = 0; i < sol.slots[t].spectrums[s].channels[c].connections.size(); i++) {
                    fprintf(file, " %d", sol.slots[t].spectrums[s].channels[c].connections[i].id);
                }
                fprintf(file, "\n");
            }
        }
    }

    fclose(file);
}

double computeConnectionThroughput(Connection &conn, int bandwidth, bool force = false) {
    if (bandwidth == 0) {
        return 0.0;
    }

    int mcs = -1;
    int maxDataRate = 11;

    if (approximatelyEqual(conn.interference, 0.0, EPS)) {
        mcs = maxDataRate;
        conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
    } else {
        double conn_SINR = (powerSender / pow(conn.distanceSR, alfa)) / (conn.interference + noise);
        conn.SINR = conn_SINR;

        // printf("%lf %.10lf %lf %lf %lf %lf %lf\n", powerSender, conn.distanceSR, alfa, pow(conn.distanceSR, alfa), conn.interference, noise, conn_SINR);

        mcs = 11;
        while (mcs > 0 && conn_SINR < SINR[mcs][bwIdx(bandwidth)]) // TODO: be careful here
            mcs--;

        // while (mcs + 1 <= maxDataRate && conn_SINR > SINR[mcs + 1][bwIdx(bandwidth)])
        //     mcs++;

        if (force) {
            printf("   ->> SINR EH %.10lf, MCS %d, BW %d\n", conn_SINR, mcs, bandwidth);
        }

        conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
    }

    return conn.throughput;
}

Channel insertInChannel(Channel newChannel, int idConn) {
    Connection conn(idConn, 0.0, 0.0, distanceMatrix[idConn][idConn]);

    for (Connection &connection : newChannel.connections) {
        connection.interference += interferenceMatrix[connection.id][conn.id];
        conn.interference += interferenceMatrix[conn.id][connection.id];
    }

    newChannel.connections.emplace_back(conn);
    newChannel.throughput = 0.0;
    for (Connection &connection : newChannel.connections) {
        computeConnectionThroughput(connection, newChannel.bandwidth);
        newChannel.throughput += connection.throughput;
    }

    return newChannel;
}

Channel deleteFromChannel(const Channel &channel, int idConn) {
    Channel newChannel(channel.bandwidth);

    for (const Connection &conn : channel.connections) {
        if (conn.id != idConn) {
            newChannel.connections.emplace_back(conn);
        }
    }

    newChannel.throughput = 0.0;
    for (Connection &conn : newChannel.connections) {
        conn.interference -= interferenceMatrix[conn.id][idConn];
        computeConnectionThroughput(conn, newChannel.bandwidth);
        newChannel.throughput += conn.throughput;
    }

    return newChannel;
}

bool reinsert(Solution &sol, Connection conn, ti3 to, bool force = false) {
    Channel &new_chan = sol.slots[get<0>(to)].spectrums[get<1>(to)].channels[get<2>(to)];
    Channel copy_new_chan = insertInChannel(new_chan, conn.id);

    double newObjective = sol.throughput - new_chan.throughput + copy_new_chan.throughput;

    bool improved = false;
    if (newObjective > sol.throughput)
        improved = true;

    if ((newObjective > sol.throughput) || force) {
        new_chan = copy_new_chan;
        //
        sol.throughput = newObjective;
    }

    return improved;
}

double computeThroughput(Solution &curr, bool force = false) {
    double OF = 0.0;

    for (int t = 0; t < curr.slots.size(); t++) {
        for (int s = 0; s < curr.slots[t].spectrums.size(); s++) {
            for (int c = 0; c < curr.slots[t].spectrums[s].channels.size(); c++) {

                double &chThroughput = curr.slots[t].spectrums[s].channels[c].throughput;
                chThroughput = 0.0;
                for (Connection &conn : curr.slots[t].spectrums[s].channels[c].connections) {
                    conn.interference = 0.0;
                    conn.throughput = 0.0;
                    for (Connection &otherConn :
                         curr.slots[t].spectrums[s].channels[c].connections) {
                        conn.interference += interferenceMatrix[conn.id][otherConn.id];
                    }
                    chThroughput += computeConnectionThroughput(
                        conn, curr.slots[t].spectrums[s].channels[c].bandwidth, force);
                }
                OF += chThroughput;
            }
        }
    }

    curr.throughput = OF;
    curr.throughput_flag = true;
    return OF;
}

void print_solution(Solution &sol, bool ok = true) {
    printf("solucao com %lu time-slots e throughput %lf\n", sol.slots.size(), sol.throughput);
    for (int t = 0; t < sol.slots.size(); t++) {
        printf("TS: %d:\n{\n", t);
        for (int s = 0; s < sol.slots[t].spectrums.size(); s++) {
            printf("Spec: %d, ", s);
            for (int c = 0; c < sol.slots[t].spectrums[s].channels.size(); c++) {
                printf("CH: %d (bw: %d)\nconns:\n", c, sol.slots[t].spectrums[s].channels[c].bandwidth);
                for (const Connection &conn : sol.slots[t].spectrums[s].channels[c].connections) {
                    printf("{%d, %lf, %lf, %lf}\n", conn.id, conn.throughput, conn.interference, conn.SINR);
                }
                puts("");
            }
        }
        printf("}\n");
    }
}

void recoverSolution(int k, int i, int j, bool clean) {
    if (clean)
        inSolution[k][i][j] = false;

    if (inSolution[k][i][j]) {
        clean = true;
    }

    if (child[k][i][j][0] != -1 && child[k][i][j][1] != -1) {
        recoverSolution(k, i, child[k][i][j][0], clean);
        recoverSolution(k, i, child[k][i][j][1], clean);
    }
}

Solution reconstruct_sol(const Solution &curr) {
    Solution ret;
    ret.slots.resize(curr.slots.size());

    for (int t = 0; t < curr.slots.size(); t++) {
        ret.slots[t].spectrums.resize(curr.slots[t].spectrums.size(),
                                      Spectrum(0.0, 0.0, std::vector<Channel>()));
    }

    for (int t = 0; t < curr.slots.size(); t++) {
        for (int s = 0; s < curr.slots[t].spectrums.size(); s++) {
            for (int c = 0; c < curr.slots[t].spectrums[s].channels.size(); c++) {
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
                    ret.slots[t].spectrums[s].channels.emplace_back(
                        curr.slots[t].spectrums[s].channels[c]);
                }
            }
        }
    }

    computeThroughput(ret);

#ifdef MDVRBSP
    compute_violation(ret);
#endif
    return ret;
}

bool allChannels20MHz(const Solution &sol) {
    for (const TimeSlot &ts : sol.slots) {
        for (const Spectrum &spectrum : ts.spectrums) {
            for (const Channel &channel : spectrum.channels) {
                if (channel.bandwidth > 20) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool fix_channels(Solution &sol) {
    bool improved = false;
    // ti3 zeroChannel = make_tuple(sol.slots.size() - 1, sol.slots[0].spectrums.size() - 1, 0);
    do {
        improved = false;
        for (int t = 0; t < sol.slots.size(); t++) {
            for (int s = 0; s < sol.slots[t].spectrums.size(); s++) {
                for (int c = 0; c < sol.slots[t].spectrums[s].channels.size(); c++) {
                    // if (make_tuple(t, s, c) == zeroChannel)
                    //     continue;

                    int idx = 0;
                    while (idx < sol.slots[t].spectrums[s].channels[c].connections.size()) {
                        Connection &conn = sol.slots[t].spectrums[s].channels[c].connections[idx];
                        double aux = computeConnectionThroughput(
                            conn, sol.slots[t].spectrums[s].channels[c].bandwidth);
                        if (approximatelyEqual(aux, 0.0, EPS)) {
                            // reinsert(sol, conn, make_tuple(t, s, c), zeroChannel, true);
                            // Remover link da tupla (t, s, c) e colocar em unscheduled
                            sol.unscheduled.emplace_back(conn.id);

                            Channel &old_chan = sol.slots[t].spectrums[s].channels[c];
                            Channel copy_old_channel = deleteFromChannel(old_chan, conn.id);

                            sol.throughput =
                                sol.throughput - old_chan.throughput + copy_old_channel.throughput;
                            old_chan = copy_old_channel;
                            improved = true;
                        } else {
                            idx++;
                        }
                    }
                }
            }
        }
    } while (improved);
    computeThroughput(sol);

#ifdef MDVRBSP
    compute_violation(sol);
#endif
    return improved;
}

Solution convertTo20MhzSol(Solution exps) {
    for (TimeSlot &ts : exps.slots) {
        for (Spectrum &sp : ts.spectrums) {
            sort(sp.channels.rbegin(), sp.channels.rend());
        }
    }

    while (!allChannels20MHz(exps)) {
        for (TimeSlot &ts : exps.slots) {
            for (Spectrum &sp : ts.spectrums) {
                for (int c = 0; c < sp.channels.size(); c++) {
                    if (sp.channels[c].bandwidth <= 20)
                        continue;

                    int newBw = sp.channels[c].bandwidth / 2;
                    Channel child1(newBw), child2(newBw);

                    std::vector<Connection> connections = sp.channels[c].connections;

                    for (Connection &conn : connections) {
                        if (rng.randInt(1)) {
                            child1.connections.emplace_back(conn);
                        } else {
                            child2.connections.emplace_back(conn);
                        }
                    }

                    sp.channels[c] = child1;
                    sp.channels.insert(sp.channels.begin() + c,
                                       child2); // TODO: it was .end(). Did it right?
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
        for (int s = 0; s < ret.slots[t].spectrums.size(); s++) {
            Spectrum &spectrum = ret.slots[t].spectrums[s];
            int curr = 0;
            while (curr < spectrum.channels.size()) {
                if (curr + 1 < spectrum.channels.size() &&
                    spectrum.channels[curr].bandwidth == spectrum.channels[curr + 1].bandwidth) {
                    Channel merged = spectrum.channels[curr];
                    merged.bandwidth *= 2;

                    for (const Connection &conn : spectrum.channels[curr + 1].connections) {
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

#ifdef MDVRBSP
    compute_violation(ret);
#endif
    return ret;
}

void K_AddDrop(Solution &sol, int K) {
    K = min(K, int(sol.unscheduled.size()));
    for (int i = 0; i < K; i++) {
        int idx = rng.randInt(sol.unscheduled.size() - 1);

        Connection conn = Connection(sol.unscheduled[idx]);

        int t = rng.randInt(sol.slots.size() - 1);
        int a = rng.randInt(sol.slots[t].spectrums.size() - 1);
        int b = rng.randInt(sol.slots[t].spectrums[a].channels.size() - 1);
        ti3 channelTo = make_tuple(t, a, b);

        reinsert(sol, conn, channelTo, true);

        //Remove connection from the unscheduled vector
        swap(sol.unscheduled[idx], sol.unscheduled.back());
        sol.unscheduled.pop_back();
    }
}

// TODO: finish this code
void K_RemoveAndInserts(Solution &sol, int K) {
    int k = 0;
    while (k < K) {
        int t = rng.randInt(sol.slots.size() - 1);
        int a = rng.randInt(sol.slots[t].spectrums.size() - 1);
        int b = rng.randInt(sol.slots[t].spectrums[a].channels.size() - 1);

        if (sol.slots[t].spectrums[a].channels[b].connections.empty())
            continue;

        k++;

        Channel &ch = sol.slots[t].spectrums[a].channels[b];
        int z = rng.randInt(ch.connections.size() - 1);
        Connection conn = ch.connections[z];

        ch = deleteFromChannel(ch, conn.id);
        sol.unscheduled.emplace_back(conn.id);
    }
    K_AddDrop(sol, K);
}

double solve(int k, int i, int j) {
    inSolution[k][i][j] = true;
    double ret = chanThroughput[k][i][j];
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

void setDP(const Solution &sol, bool ok = false) {
    memset(chanThroughput, 0, sizeof chanThroughput);
    memset(inSolution, false, sizeof inSolution);

    for (int t = 0; t < sol.slots.size(); t++) {
        for (int s = 0; s < sol.slots[t].spectrums.size(); s++) {
            for (int c = 0; c < sol.slots[t].spectrums[s].channels.size(); c++) {
                chanThroughput[t][s][c] = sol.slots[t].spectrums[s].channels[c].throughput;
            }
        }
    }
}

double calcDP(const Solution &sol, bool ok = false) {
    double OF = 0.0;
    for (int t = 0; t < sol.slots.size(); t++) {
        for (int s = 0; s < sol.slots[t].spectrums.size(); s++) {
            for (int c = 0; c < sol.slots[t].spectrums[s].channels.size(); c++) {
                if (parent[t][s][c] == -1) {
                    double ret = solve(t, s, c);
                    OF += ret;
                }
            }
        }
    }

    return OF;
}

bool is_okay(Solution &sol) {
    set<int> saw;
    for (TimeSlot &ts : sol.slots) {
        for (Spectrum &sp : ts.spectrums) {
            for (Channel &ch : sp.channels) {
                for (Connection &conn : ch.connections) {
                    if (saw.count(conn.id)) {
                        printf("ja vi %d\n", conn.id);
                        return false;
                    } else {
                        saw.insert(conn.id);
                    }
                }
            }
        }
    }

    return true;
}

bool is_okay2(Solution &sol) {
    for (TimeSlot &ts : sol.slots) {
        for (Spectrum &sp : ts.spectrums) {
            for (Channel &ch : sp.channels) {
                set<int> saw;
                for (Connection &conn : ch.connections) {
                    if (saw.count(conn.id)) {
                        printf("ja vi %d\n", conn.id);
                        print_solution(sol);
                        return false;
                    } else {
                        saw.insert(conn.id);
                    }
                }
            }
        }
    }

    return true;
}

Solution local_search(Solution &multiple, Solution &curr) {
    bool improved = false;

    do {
        improved = false;
        for (int i = 0; i < n_connections; i++) {
            Solution mult_clean(multiple);
            for (int t = 0; t < mult_clean.slots.size(); t++) {
                for (Spectrum &sp : mult_clean.slots[t].spectrums) {
                    for (Channel &channel : sp.channels) {
                        for (const Connection &conn : channel.connections) {
                            if (conn.id == i) {
                                channel = deleteFromChannel(channel, conn.id);
                                break;
                            }
                        }
                    }
                }
            }

            Solution mult_cont(mult_clean);
            for (int t = 0; t < mult_cont.slots.size(); t++) {
                for (Spectrum &sp : mult_cont.slots[t].spectrums) {
                    for (Channel &ch : sp.channels) {
                        ch = insertInChannel(ch, i);
                    }
                }
            }

            double bestOF = -1;
            ti3 best_ch = make_tuple(-1, -1, -1);

            for (int t = 0; t < curr.slots.size(); t++) {
                for (int s = 0; s < curr.slots[t].spectrums.size(); s++) {
                    for (int c = 0; c < curr.slots[t].spectrums[s].channels.size(); c++) {
                        setDP(mult_clean);
                        int c_ch = c;
                        while (c_ch != -1) {
                            chanThroughput[t][s][c_ch] =
                                mult_cont.slots[t].spectrums[s].channels[c_ch].throughput;
                            c_ch = parent[t][s][c_ch];
                        }

                        double OF = calcDP(multiple);
                        if (OF > bestOF) {
                            bestOF = OF;
                            best_ch = make_tuple(t, s, c);
                        }
                    }
                }
            }

            if (bestOF > curr.throughput) {
                // printf("entrei %lf %lf\n", bestOF, curr.throughput);
                improved = true;

                curr.throughput = bestOF;
                for (int t = 0; t < multiple.slots.size(); t++) {
                    for (int s = 0; s < multiple.slots[t].spectrums.size(); s++) {
                        for (int c = 0; c < multiple.slots[t].spectrums[s].channels.size(); c++) {
                            Channel &ch = multiple.slots[t].spectrums[s].channels[c];
                            for (const Connection &conn : ch.connections) {
                                if (conn.id == i) {
                                    ch = deleteFromChannel(ch, conn.id);
                                    break;
                                }
                            }
                        }
                    }

                    int new_t = get<0>(best_ch);
                    int new_sp = get<1>(best_ch);
                    int new_ch = get<2>(best_ch);
                    while (new_ch != -1) {
                        multiple.slots[new_t].spectrums[new_sp].channels[new_ch] = insertInChannel(
                            multiple.slots[new_t].spectrums[new_sp].channels[new_ch], i);
                        new_ch = parent[new_t][new_sp][new_ch];
                    }
                }
            }
        }

    } while (improved);

    Solution to_ret = reconstruct_sol(multiple);
    return to_ret;
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

void rawInsert(Solution &sol, int conn, ti3 to) {
    sol.slots[get<0>(to)].spectrums[get<1>(to)].channels[get<2>(to)].connections.emplace_back(
        Connection(conn));
}

vector<Channel> split(Channel toSplit) {
    int newBw = toSplit.bandwidth / 2;
    vector<Channel> ret;
    ret.emplace_back(Channel(0.0, 0.0, newBw, vector<Connection>()));
    ret.emplace_back(Channel(0.0, 0.0, newBw, vector<Connection>()));

    vector<Connection> connToInsert = toSplit.connections;
    sort(connToInsert.rbegin(), connToInsert.rend());

    double bestThroughputSoFar = 0.0;
    for (int i = 0; i < connToInsert.size(); i++) {
        double bestThroughputIteration = 0.0, testando = 0.0;
        Channel bestChannel(0.0, 0.0, newBw, vector<Connection>());
        int bestIdx = -1;

        for (int c = 0; c < ret.size(); c++) {
            Channel inserted = insertInChannel(ret[c], connToInsert[i].id);
            double resultThroughput = inserted.throughput;

            double g_throughput = (bestThroughputSoFar - ret[c].throughput) + resultThroughput;
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
        bestThroughputSoFar = bestThroughputIteration;
    }
    computeChannelsThroughput(ret);
    assert(approximatelyEqual(ret[0].throughput + ret[1].throughput, bestThroughputSoFar));
    return ret;
}

void insertConnectionAt(int ts, int sp, int ch, int conn_id, Solution &curr) {}

void remove_timeslot(Solution &curr) {
    int ts_idx = -1;
    double ts_inter = -1.0;
    for (int t = 0; t < int(curr.slots.size()); t++) {
        if (curr.slots[t].interference > ts_inter) {
            ts_idx = t;
            ts_inter = curr.slots[t].interference;
        }
    }

    // Removo o time_slot
    TimeSlot cp_ts(curr.slots[ts_idx]);
    swap(curr.slots[ts_idx], curr.slots.back());
    curr.slots.pop_back();

    // Redestribuir conexoes
    for (const Spectrum &sp : cp_ts.spectrums) {
        for (const Channel &ch : sp.channels) {
            for (const Connection &conn : ch.connections) {
                int n_ts = rng.randInt(curr.slots.size());
                int n_sp = rng.randInt(curr.slots[n_ts].spectrums.size());
                int n_ch = rng.randInt(curr.slots[n_ts].spectrums[n_ch].channels.size());
                // inserir conexao no local sorteado e re-computar throughput de todo mundo
                insertConnectionAt(n_ts, n_sp, n_ch, conn.id, curr);
            }
        }
    }
}

inline void perturbation(Solution &sol, int kkmul) {
    int rnd = rng.randInt(1);
    if (rnd) {
        // puts("entrei 1");
        K_AddDrop(sol, kkmul);
    } else {
        // puts("entrei 2");
        K_RemoveAndInserts(sol, kkmul);
    }

    fix_channels(sol);
}

void init(void);

Solution vns(void);

Solution constructive_heuristic(void);

bool is_feasible(const Solution &);

#endif
