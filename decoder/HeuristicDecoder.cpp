//
// Created by José Joaquim on 03/03/20.
//

#include "HeuristicDecoder.h"

using namespace std;

const int X_c = 0;
const int Y_c = 1;
const double EPS = 1e-7;

int nConnections, nSpectrums;
vector<vector<double>> dataRates, SINR, beta;
double distanceMatrix[MAX_CONN][MAX_CONN], interferenceMatrix[MAX_CONN][MAX_CONN];
double senders[MAX_CONN][2], receivers[MAX_CONN][2];
static double affectance[MAX_CONN][MAX_CONN];
double powerSender, alfa, noise, ttm;

static vector<Spectrum> init_conf;
static vector<int> spectrum_size;

random_device rd;
pair<int, int> zeroChannel;

auto whatever = std::default_random_engine{rd()};
MTRand rng;

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

inline double distance(double X_si, double Y_si, double X_ri, double Y_ri) {
    return hypot((X_si - X_ri), (Y_si - Y_ri));
}

void distanceAndInterference() {
    for (int i = 0; i < nConnections; i++) {
        double X_si = receivers[i][X_c];
        double Y_si = receivers[i][Y_c];

        for (int j = 0; j < nConnections; j++) {

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

    // for (int i = 0; i < nConnections; i++) {
    //     for (int j = 0; j < nConnections; j++) {
    //         printf("%.10lf ", interferenceMatrix[i][j]);
    //     }
    //     puts("");
    // }
}

double convertDBMToMW(double _value) {
    double result = 0.0, b;

    b = _value / 10.0;     // dBm dividido por 10
    result = pow(10.0, b); // Converte de DBm para mW

    return result;
}

void convertTableToMW(const vector<vector<double>> &_SINR, vector<vector<double>> &_SINR_Mw) {
    double result, b;
    for (int i = 0; i < _SINR_Mw.size(); i++) {
        for (int j = 0; j < _SINR_Mw[i].size(); j++) {

            if (_SINR[i][j] != 0) {
                b = _SINR[i][j] / 10.0; // dBm divided by 10
                result = pow(10.0, b);  // Convert DBM to mW

                _SINR_Mw[i][j] = result;
            } else {
                _SINR_Mw[i][j] = 0;
            }
        }
    }
}

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

    zeroChannel = {init_conf.size(), 0};
}

void loadData() {
    scanf("%d %lf %lf %lf %d", &nConnections, &alfa, &noise, &powerSender, &nSpectrums);

    for (int i = 0; i < nSpectrums; i++) {
        int s;
        scanf("%d", &s);
        spectrum_size.emplace_back(s);
        init_conf.emplace_back(s, 0, vector<Channel>());
    }

    if (noise != 0) {
        noise = convertDBMToMW(noise);
    }

    for (int i = 0; i < nConnections; i++) {
        double a, b;
        scanf("%lf %lf", &a, &b);
        receivers[i][0] = a;
        receivers[i][1] = b;
    }

    for (int i = 0; i < nConnections; i++) {
        double a, b;
        scanf("%lf %lf", &a, &b);
        senders[i][0] = a;
        senders[i][1] = b;
    }

    for (int i = 0; i < nConnections; i++) {
        double bt;
        scanf("%lf", &bt);
#ifdef MDVRBSP
        gma.emplace_back(bt);
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

    for (int i = 0; i < nConnections; i++) {
        for (int j = 0; j < nConnections; j++) {
            double value = powerSender / pow(distanceMatrix[i][j], alfa);
            affectance[i][j] = value;
        }
    }

#ifdef MDVRBSP
    // Select random dataRates from the table
    vector<double> aux;
    for (int i = 0; i < int(dataRates.size()); i++)
        for_each(dataRates[i].begin(), dataRates[i].end(),
                 [&aux](double x) { aux.emplace_back(x); });

    for (int i = 0; i < nConnections; i++)
        gma[i] = aux[rng.randInt(aux.size() - 1)];

    beta.assign(nConnections, vector<double>(4, 0));
    for (int i = 0; i < nConnections; i++) {

        for (int j = 0; j < 4; j++) {
            double value = gammaToBeta(gma[i], j);
            beta[i][j] = value;
        }
    }
#endif
}

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

bool double_equals(double a, double b, double epsilon) { return std::abs(a - b) < epsilon; }

bool operator>(const Solution &o1, const Solution &o2) { return operator<(o2, o1); }

bool operator<(const Solution &o1, const Solution &o2) {
    assert(o1.throughputFlag && o2.throughputFlag);
    return o1.totalThroughput < o2.totalThroughput;
}

double computeConnectionThroughput(Connection &conn, int bandwidth, bool force) {
    if (bandwidth == 0)
        return 0.0;

    int mcs = -1;
    int maxDataRate = bandwidth == 20 ? 8 : 9;

    if (double_equals(conn.interference, 0.0)) {
        mcs = maxDataRate;
        conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
    } else {
        double conn_SINR = (powerSender / pow(conn.distanceSR, alfa)) / (conn.interference + noise);
        conn.SINR = conn_SINR;

        while (mcs + 1 <= maxDataRate && conn_SINR > SINR[mcs + 1][bwIdx(bandwidth)])
            mcs++;

        if (force) {
            printf("   ->> SINR EH %.10lf, MCS %d, BW %d\n", conn_SINR, mcs, bandwidth);
        }

        conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
    }

    return conn.throughput;
}

void testConnection(Connection &conn, int bandwidth, bool force = true) {
    printf("conn %d:\n", conn.id);
    printf("  interference: %.10lf\n", conn.interference);
    printf("  bw: %d\n", bandwidth);
    printf("  distanceSR: %.10lf\n", conn.distanceSR);
    computeConnectionThroughput(conn, bandwidth, true);
}

Channel insertInChannel(Channel newChannel, int idConn) {
    Connection conn(idConn, 0.0, 0.0, distanceMatrix[idConn][idConn]);

    for (Connection &connection : newChannel.connections) {
        //    connection.interference += interferenceMatrix[conn.id][connection.id];
        //    conn.interference += interferenceMatrix[connection.id][conn.id];
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

double computeThroughput(Solution &curr, bool force) {
    double OF = 0.0;

    for (int s = 0; s < curr.spectrums.size(); s++) {
        for (int c = 0; c < curr.spectrums[s].channels.size(); c++) {
            if (make_pair(s, c) == zeroChannel)
                break;

            double &chThroughput = curr.spectrums[s].channels[c].throughput;
            chThroughput = 0.0;
            for (Connection &conn : curr.spectrums[s].channels[c].connections) {
                conn.interference = 0.0;
                conn.throughput = 0.0;
                for (Connection &otherConn : curr.spectrums[s].channels[c].connections) {
                    conn.interference += interferenceMatrix[conn.id][otherConn.id];
                }
                chThroughput += computeConnectionThroughput(
                    conn, curr.spectrums[s].channels[c].bandwidth, force);
            }
            OF += chThroughput;
        }
    }

    curr.totalThroughput = OF;
    curr.throughputFlag = true;
    return OF;
}

void insertInSpectrum(Solution &sol, vector<Channel> &channels, int specId) {
    sol.throughputFlag = false;
    for (const Channel &ch : channels) {
        sol.spectrums[specId].channels.emplace_back(ch);
    }

    computeThroughput(sol);
}

void rawInsert(Solution &sol, int conn, ii where) {
    sol.spectrums[where.first].channels[where.second].connections.emplace_back(Connection(conn));
}

void rawRemove(Solution &sol, int conn, ii where) {
    for (int idx = 0; idx < sol.spectrums[where.first].channels[where.second].connections.size();
         idx++) {
        if (sol.spectrums[where.first].channels[where.second].connections[idx].id == conn) {
            sol.spectrums[where.first].channels[where.second].connections.erase(
                sol.spectrums[where.first].channels[where.second].connections.begin() + idx);
            return;
        }
    }
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

double some(const Connection &c1, const Connection &c2) { return c1.throughput + c2.throughput; }

vector<Channel> split(Channel toSplit, int spectrum) {
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
            bool two = double_equals(g_throughput, bestThroughputIteration) &&
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
    assert(double_equals(ret[0].throughput + ret[1].throughput, bestThroughputSoFar));
    return ret;
}

void printChannel(const Channel &chan, ii where = {-1, -1}) {
    printf("Atual estado do canal {%d, %d} (%d MHz):\n", where.first, where.second, chan.bandwidth);
    printf("   Throughput: %.4lf\n", chan.throughput);
    printf("   Conexoes:");
    double sum = 0.0;
    for (const Connection &conn : chan.connections) {
        sum += conn.throughput;
        printf(" {%d, %.4lf, %.4lf} ", conn.id, conn.distanceSR, conn.throughput);
    }
    puts("");

    assert(double_equals(sum, chan.throughput));
}

Connection::Connection(int id, double throughput, double interference, double distanceSR)
    : id(id), throughput(throughput), interference(interference), distanceSR(distanceSR) {}

Channel::Channel(double throughput, double interference, int bandwidth,
                 const vector<Connection> &connections)
    : throughput(throughput), interference(interference), bandwidth(bandwidth),
      connections(connections) {}

Channel::Channel(int bandwidth) : bandwidth(bandwidth) {
    interference = 0.0;
    throughput = 0.0;
    connections = vector<Connection>();
}

Spectrum::Spectrum(int maxFrequency, int usedFrequency, const vector<Channel> &channels)
    : maxFrequency(maxFrequency), usedFrequency(usedFrequency), channels(channels) {}

Solution::Solution(const vector<Spectrum> &spectrums, double total, bool flag)
    : spectrums(spectrums), totalThroughput(total), throughputFlag(flag) {}

Solution::Solution() {
    throughputFlag = true;
    totalThroughput = 0.0;
}

Connection::Connection(int id) : id(id) {
    interference = 0.0;
    throughput = 0.0;
    distanceSR = distanceMatrix[id][id];
}

bool addChannel(Solution &sol, int spec, int bw) {
    int freeFrequency = sol.spectrums[spec].maxFrequency - sol.spectrums[spec].usedFrequency;
    if (freeFrequency < bw) {
        return false;
    }

    Channel newChannel(bw);
    sol.spectrums[spec].channels.emplace_back(newChannel);
    sol.spectrums[spec].usedFrequency += bw;
    return true;
}

bool removeChannel(Solution &sol, ii where) {
    double usedBw = sol.spectrums[where.first].channels[where.second].bandwidth;
    sol.spectrums[where.first].usedFrequency -= usedBw;
    sol.spectrums[where.first].channels.erase(sol.spectrums[where.first].channels.begin() +
                                              where.second);
    return true;
}

void Solution::printSolution(FILE *file) {
    if (file == nullptr)
        file = stdout;

    int cont = 0;
    int arr[] = {24, 36, 42, 44};
    // TODO
    // for (int s = 0; s < spectrums.size() - 1; s++) {
    //     const Spectrum &sc = spectrums[s];
    //     for (const Channel &ch : sc.channels) {
    //         for (const Connection &conn : ch.connections) {
    //             cont++;
    //             fprintf(file, "%d %d %d %d %lf\n", conn.id, arr[bwIdx(ch.bandwidth)],
    //                     bwIdx(ch.bandwidth), computeConnectionMCS(conn, ch.bandwidth),
    //                     conn.interference);
    //         }
    //         arr[bwIdx(ch.bandwidth)]--;
    //     }
    // }

    fprintf(file, "TOTAL OF %d CONNECTIONS\n", cont);
}

int convertChromoBandwidth(double value) {
    if (value < .25)
        return 20;
    else if (value < .5)
        return 40;
    else if (value < .75)
        return 80;

    return 160;
}

int insertFreeChannel(Solution &sol, int conn, int band, vector<double> &variables) {
    int freeSpec = 0, freeSpec2 = 0, idxSpectrum = -1;
    int auxBw = -1;
    bool inserted = false;
    for (int s = 0; s < sol.spectrums.size(); s++) {
        freeSpec = sol.spectrums[s].maxFrequency - sol.spectrums[s].usedFrequency;

        if (freeSpec >= band) {
            sol.spectrums[s].usedFrequency += band;

            Channel newChannel = insertInChannel(Channel(band), conn);
            sol.spectrums[s].channels.emplace_back(newChannel);

            sol.totalThroughput += newChannel.throughput;
            inserted = true;
            return band;
        } else if (freeSpec > freeSpec2) { // TODO: what is the condition?
            freeSpec2 = freeSpec;
            idxSpectrum = s;

            if (freeSpec2 >= 80) {
                auxBw = 80;
            } else if (freeSpec2 >= 40) {
                auxBw = 40;
            } else {
                auxBw = 20;
            }
        }
    }

    assert(auxBw != -1);
    if (!inserted && freeSpec2 > 0) {
        Channel newChannel = insertInChannel(Channel(auxBw), conn);
        sol.spectrums[idxSpectrum].usedFrequency += auxBw;
        sol.spectrums[idxSpectrum].channels.emplace_back(newChannel);
        sol.totalThroughput += newChannel.throughput;

        if (auxBw != band) {
            switch (newChannel.bandwidth) {
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

    return auxBw;
}

void insertBestChannel(Solution &sol, int conn, int band, vector<double> &variables) {
    // TODO: remind to check the output of this function.
    double currentThroughput = sol.totalThroughput, bestThroughputIteration = sol.totalThroughput;
    bool inserted = false;
    Channel newChannel(band);
    ii nCh = {-1, -1};
    for (int s = 0; s < sol.spectrums.size(); s++) {
        for (int c = 0; c < sol.spectrums[s].channels.size(); c++) {
            Channel channelInsert = insertInChannel(sol.spectrums[s].channels[c], conn);

            double auxThroughput = currentThroughput - sol.spectrums[s].channels[c].throughput +
                                   channelInsert.throughput;

            if (auxThroughput > bestThroughputIteration) {
                bestThroughputIteration = auxThroughput;
                newChannel = channelInsert;
                inserted = true;
                nCh = {s, c};
            }
        }
    }

    if (inserted) {
        sol.totalThroughput = sol.totalThroughput -
                              sol.spectrums[nCh.first].channels[nCh.second].throughput +
                              newChannel.throughput;
        sol.spectrums[nCh.first].channels[nCh.second] = newChannel;

        if (newChannel.bandwidth != band) {
            switch (newChannel.bandwidth) {
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
}

double buildVRBSPSolution(vector<double> variables, vector<int> permutation) {
    int totalSpectrum = 160 + 240 + 100, totalUsedSpectrum = 0.0;

    vector<Channel> auxCh;
    Spectrum spec1(160, 0, auxCh);
    Spectrum spec2(240, 0, auxCh);
    Spectrum spec3(100, 0, auxCh);
    Solution sol({spec1, spec2, spec3}, 0.0, true);

    // First, insert in free channels
    int idx = 0;
    while (idx < permutation.size() && totalUsedSpectrum < totalSpectrum) {
        int connection = permutation[idx] / 2;
        int bandWidth = convertChromoBandwidth(variables[permutation[idx] + 1]);
        totalUsedSpectrum += insertFreeChannel(sol, connection, bandWidth, variables);
        idx++;
    }

    // Second, insert in the best channels
    while (idx < permutation.size()) {
        int connection = permutation[idx] / 2;
        int bandWidth = convertChromoBandwidth(variables[permutation[idx] + 1]);
        insertBestChannel(sol, connection, bandWidth, variables);
        idx++;
    }

    double throughput = sol.totalThroughput;
    return -1.0 * throughput;
}

double Solution::decode(std::vector<double> variables) const {
    double fitness = 0.0;

    vector<pair<double, int>> ranking;
    for (int i = 0; i < variables.size(); i += 2)
        ranking.emplace_back(variables[i], i);

    sort(ranking.begin(), ranking.end());

    vector<int> permutation;
    for (int i = 0; i < ranking.size(); i++)
        permutation.emplace_back(ranking[i].second);

    fitness = buildVRBSPSolution(variables, permutation);
    return fitness;
}

double computeConnectionThroughput(Connection &conn, int bandwidth) {
    if (bandwidth == 0)
        return 0.0;

    int mcs = -1;
    int maxDataRate = 11;

    if (approximatelyEqual(conn.interference, 0.0, EPS)) {
        mcs = maxDataRate;
        conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
    } else {
        double conn_SINR = (powerSender / pow(conn.distanceSR, alfa)) / (conn.interference + noise);
        conn.SINR = conn_SINR;

        mcs = 11;
        while (mcs >= 0 && conn_SINR < SINR[mcs][bwIdx(bandwidth)])
            mcs--;

        if (mcs < 0)
            conn.throughput = 0.0;
        else
            conn.throughput = dataRates[mcs][bwIdx(bandwidth)];
    }

    return conn.throughput;
}
