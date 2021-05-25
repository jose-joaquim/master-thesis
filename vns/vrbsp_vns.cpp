#include "VNS.h"

using namespace std;

Solution constructive_heuristic() {
    Solution ret(init_conf, 0.0, true);
    vector<int> links;
    for (int i = 0; i < n_connections; i++)
        links.emplace_back(i);

    shuffle(links.begin(), links.end(), whatever);

    Solution retCopy(ret);
    vector<int> zeroConn;
    double bestThroughputSoFar = 0.0;
    for (int conn : links) {
        double bestThroughputIteration = bestThroughputSoFar;
        vector<Channel> bestSplitChannels;
        Channel bestChannel(0.0, 0.0, 0, vector<Connection>());
        ii where = {-1, -1};
        bool isSplit = false, inserted = false;
        int bestBandwidthSplit = -1, bandwidthSplit = -1;

        for (int s = 0; s < retCopy.slots[0].spectrums.size(); s++) {
            for (int c = 0; c < retCopy.slots[0].spectrums[s].channels.size(); c++) {
                const Channel &currentChannel = retCopy.slots[0].spectrums[s].channels[c];
                bandwidthSplit = -1;

                Channel channelInserted = insertInChannel(currentChannel, conn);
                vector<Channel> channelSplit;
                if (channelInserted.bandwidth >= 40) {
                    channelSplit = split(currentChannel);

                    Channel split0 = insertInChannel(channelSplit[0], conn);
                    Channel split1 = insertInChannel(channelSplit[1], conn);

                    int op = 0;
                    double conf1 = split0.throughput + channelSplit[1].throughput;
                    double conf2 = split1.throughput + channelSplit[0].throughput;
                    if (conf1 > conf2) {
                        channelSplit[0] = split0;
                    } else if (conf2 > conf1) {
                        op = 1;
                        channelSplit[1] = split1;
                    } else if (split0.connections.size() < split1.connections.size()) {
                        op = 2;
                        channelSplit[0] = split0;
                    } else {
                        op = 3;
                        channelSplit[1] = split1;
                    }
                }

                double ammount = channelInserted.throughput;
                bool betterSplit = false;
                if (!channelSplit.empty()) {
                    ammount = max(channelSplit[0].throughput + channelSplit[1].throughput, ammount);
                    betterSplit = (channelSplit[0].throughput + channelSplit[1].throughput) >
                                  channelInserted.throughput;
                }

                bandwidthSplit =
                    betterSplit ? currentChannel.bandwidth / 2 : currentChannel.bandwidth;
                double totalThroughputIteration =
                    bestThroughputSoFar - currentChannel.throughput + ammount;

                bool cond1 = totalThroughputIteration > bestThroughputIteration;
                bool cond2 =
                    approximatelyEqual(totalThroughputIteration, bestThroughputIteration) &&
                    bandwidthSplit > bestBandwidthSplit;
                bool cond3 =
                    approximatelyEqual(totalThroughputIteration, bestThroughputIteration) &&
                    (bandwidthSplit == bestBandwidthSplit) && betterSplit == false &&
                    isSplit == true;

                if (cond1 || cond2 || cond3) {
                    bestThroughputIteration = totalThroughputIteration;
                    inserted = true;
                    where = {s, c};

                    bestBandwidthSplit = bandwidthSplit;
                    if (betterSplit) {
                        isSplit = true;
                        bestSplitChannels = channelSplit;
                    } else {
                        isSplit = false;
                        bestChannel = channelInserted;
                    }
                }
            }
        }

        double aux = bestThroughputSoFar;
        bestThroughputSoFar = bestThroughputIteration;

        if (inserted) {
            if (isSplit) {
                retCopy.slots[0].spectrums[where.first].channels.erase(
                    retCopy.slots[0].spectrums[where.first].channels.begin() + where.second);
                retCopy.slots[0].spectrums[where.first].channels.insert(
                    retCopy.slots[0].spectrums[where.first].channels.begin() + where.second,
                    bestSplitChannels.begin(), bestSplitChannels.end());
            } else {
                retCopy.slots[0].spectrums[where.first].channels[where.second] = bestChannel;
            }
            computeThroughput(retCopy);

            if (!approximatelyEqual(retCopy.throughput, bestThroughputSoFar)) {
                printf("vixe %.3lf %.3lf\n", retCopy.throughput, bestThroughputSoFar);
            }
            assert(approximatelyEqual(retCopy.throughput, bestThroughputSoFar));
        } else {
            zeroConn.push_back(conn);
        }
    }

    computeThroughput(retCopy);
    assert(approximatelyEqual(bestThroughputSoFar, retCopy.throughput));

    for (int i = 0; i < zeroConn.size(); i++) {
        retCopy.unscheduled.emplace_back(zeroConn[i]);
    }

    ret = retCopy;
    return ret;
}

Solution vns(string filePrefix) {
    // 1st step: constructive heuristic
    // 2nd step: dp algorithm
    // 3rd step: local search
    // 4th step: main loop until stopping condition
    //         : 1) perturbation
    //         : 2) local search
    //         : 3) update

    // ~~~~~~~~~~~~~~~~~~~~~~~`
    // TODO: ===============>> O THROUGHPUT DE UMA SOLUCAO EH O RETORNADO PELA DP!!!
    // ~~~~~~~~~~~~~~~~~~~~~~~`

    FILE *file_comparative = nullptr;
    if (!filePrefix.empty())
        file_comparative = fopen(filePrefix.c_str(), "w");

    Solution init_sol = constructive_heuristic(); // TODO (?)
    Solution delta = convertTo20MhzSol(init_sol); // DONE

    Solution rep = multipleRepresentation(delta); // DONE
    setDP(rep);
    double retOF = calcDP(rep);

    // Then, reconstruct optimal local solution
    Solution incumbent = reconstruct_sol(rep); // DONE
    incumbent.throughput = init_sol.throughput;

    Solution local_max = incumbent;
    double old_value = incumbent.throughput;
    fprintf(file_comparative, "%lf %lf\n", 0.0, old_value);

    int K_MUL = max(1, n_connections / 100);
    int K_MAX = 10;
    startTime = clock();
    while (!stop()) {
        int k = 1;
        while (!stop() && k <= K_MAX) {
            delta = local_max;

            perturbation(delta, k * K_MUL);                    // DONE
            Solution multiple = multipleRepresentation(delta); // DONE

            setDP(multiple);
            delta.throughput = calcDP(multiple);

            // printf("%lf %lf\n", incumbent.throughput, delta.throughput);

            Solution explicit_sol = local_search(multiple, delta); // TODO
            fix_channels(explicit_sol);                            // DONE

            delta = convertTo20MhzSol(explicit_sol); // DONE
            assert(essentiallyEqual(explicit_sol.throughput, delta.throughput));
            // printf("opa %lf %lf\n", delta.throughput, local_max.throughput);
            if (definitelyGreaterThan(delta.throughput, local_max.throughput)) {
                // printf("delta better %lf => %lf\n", local_max.throughput, delta.throughput);
                k = 1;
                local_max = delta;
            } else {
                k += 1;
            }

            if (definitelyGreaterThan(local_max.throughput, incumbent.throughput)) {
                double elapsed_time = (((double)(clock() - startTime)) / CLOCKS_PER_SEC);
                if (file_comparative != nullptr)
                    fprintf(file_comparative, "%lf %lf\n", elapsed_time, local_max.throughput);

                // printf("melhorei! %lf %lf => %lf\n", elapsed_time, incumbent.throughput,
                // local_max.throughput);
                incumbent = explicit_sol;
            }
        }
    }

    // double output_value = incumbent.throughput;
    // printf("%.3lf %.3lf\n", old_value, output_value);
    return incumbent;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        puts("argument error");
        return 1;
    }

    string path_input = "../instances/md-vrbsp/U_";
    path_input += string(argv[1]);
    path_input += "/MD-VRBSP_U_";
    path_input += string(argv[1]);
    path_input += "_";
    path_input += string(argv[2]);
    path_input += ".txt";

    freopen(path_input.c_str(), "r", stdin);
    maximumTime = stoi(argv[4]) * 1.0;
    read_data(); // Input will be redirected in the command call

    string objective_improvements = string(argv[3]);
    objective_improvements += "/objective_improvements" + string(argv[2]);
    objective_improvements += ".txt";

    Solution aux = vns(objective_improvements);

    // print_solution(aux);
    string path_out = string(argv[3]);
    path_out += "/solution" + string(argv[2]);
    path_out += ".txt";
    cout << path_out << endl;
    FILE *file_out = fopen(path_out.c_str(), "w");
    print_solution_to_file(aux, file_out);

    string path_obj = string(argv[3]);
    path_obj += "/objectives.txt";
    cout << path_obj << endl;
    FILE *file_obj = fopen(path_obj.c_str(), "a");
    print_objective_to_file(aux, file_obj);
    return 0;
}
