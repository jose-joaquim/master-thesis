// run 8 1 results/mdvrbsp/U_8 10
#include "VNS.h"

bool is_feasible(const Solution &ret, bool retOption = true) {
    int cnt = 0;
    set<int> scheduled;
    for (const TimeSlot &ts : ret.slots) {
        for (const Spectrum &sp : ts.spectrums) {
            for (const Channel &ch : sp.channels) {
                for (const Connection &conn : ch.connections) {
                    // printf("conn %d\n", conn.id);
                    cnt += 1;
                    scheduled.insert(conn.id);
                    if (conn.throughput < gma[conn.id]) {
                        printf("ops... conn %d %lf %lf\n", conn.id, conn.throughput, gma[conn.id]);
                        return false;
                    }
                }
            }
        }
    }

    if (cnt != n_connections) {
        for (auto x : scheduled)
            printf("%d ", x);

        puts("");
    }

    if (retOption) {
        assert(cnt == n_connections);
        return cnt == n_connections;
    }

    assert(cnt >= n_connections);
    return cnt >= n_connections;
}

bool can_split(const Channel &ch) {
    if (ch.bandwidth <= 40)
        return false;

    int new_bw = ch.bandwidth / 2;
    for (const Connection &conn : ch.connections) {
        if (gma[conn.id] > dataRates[11][bwIdx(new_bw)]) {
            return false;
        }
    }

    return true;
}

bool try_insert(int conn_id, Channel &ch) {
    ch = insertInChannel(ch, conn_id);

    for (Connection &conn : ch.connections)
        if (definitelyLessThan(conn.throughput, gma[conn.id]))
            return false;

    return true;
}

Solution constructive_heuristic(string filePrefix) {
    TimeSlot dummy_ts(init_conf);
    Solution ret({dummy_ts});

    vector<int> links;
    for (int i = 0; i < n_connections; i++)
        links.emplace_back(i);

    shuffle(links.begin(), links.end(), whatever);
    for (int conn : links) {
        bool success = false;
        for (int t = 0; t < ret.slots.size(); t++) {
            for (int s = 0; s < ret.slots[t].spectrums.size(); s++) {
                if (success)
                    break;

                for (int c = 0; c < ret.slots[t].spectrums[s].channels.size(); c++) {
                    if (make_tuple(t, s, c) == zeroChannel)
                        continue;

                    // 1. Tentar inserir no canal c
                    //    - sucesso: break
                    //    - falha: verificar se eh possivel split(c)
                    //       - se sim, tentar inserir em split(c)
                    //       - sucesso: break
                    if (success)
                        break;

                    Channel cp_channel = ret.slots[t].spectrums[s].channels[c];
                    success = try_insert(conn, cp_channel);

                    if (success) {
                        // printf(" consegui inserir %d sem quebrar o canal\n", conn);
                        ret.slots[t].spectrums[s].channels[c] = cp_channel;
                        break;
                    }

                    if (can_split(ret.slots[t].spectrums[s].channels[c])) {
                        vector<Channel> ch_split =
                            split(ret.slots[t].spectrums[s].channels[c]); // TODO

                        cp_channel = ch_split[0];
                        success = try_insert(conn, cp_channel);

                        if (success) {
                            // puts("opa 1");
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
                            // puts("opa 2");
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

            // if (!success) {
            //     TimeSlot new_ts(dummy_ts);
            // }

            // 2. Se chegou ate aqui e link nao foi inserido:
            //    - crie novo time-slot e insira no canal de maior largura de banda
            if (!success) {
                TimeSlot new_ts(dummy_ts);
                for (int s = 0; s < new_ts.spectrums.size() && !success; s++) {
                    for (int c = 0; c < new_ts.spectrums[s].channels.size() && !success; c++) {
                        int bw = new_ts.spectrums[s].channels[c].bandwidth;
                        if (bw >= 160) {
                            Connection new_conn(conn);
                            new_conn.SINR = 1000000007;
                            new_conn.throughput = dataRates[11][bwIdx(bw)];
                            new_ts.spectrums[s].channels[c].connections.push_back(new_conn);
                            ret.slots.emplace_back(new_ts);
                            success = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    vector<Channel> aux;
    aux.emplace_back(Channel());
    ret.slots[0].spectrums.emplace_back(0, 0, aux);

    computeThroughput(ret);

    if (!filePrefix.empty()) {
        FILE *objFile = fopen(filePrefix.c_str(), "w");
        fprintf(objFile, "%.3lf %lu\n", (clock() - startTime) / CLOCKS_PER_SEC, ret.slots.size());
    }
    assert(is_feasible(ret));
    return ret;
}

Solution delete_time_slot(const Solution &sol) {
    Solution ret = Solution(sol);

    assert(ret.slots[0].spectrums[3].channels[0].connections.empty());
    int to_del = rng.randInt(ret.slots.size() - 1);
    vector<Connection> links;
    for (const Spectrum &sp : ret.slots[to_del].spectrums)
        for (const Channel &ch : sp.channels)
            for (const Connection &conn : ch.connections)
                links.emplace_back(conn);

    swap(ret.slots[to_del], ret.slots.back());
    ret.slots.pop_back();

    if (to_del == 0)
        ret.slots[0].spectrums.emplace_back(0, 0, vector<Channel>{Channel(0, 0, 0, 0, links)});
    else
        ret.slots[0].spectrums[3].channels[0].connections = links;

    K_AddDrop(ret, links.size());
    computeViolation(ret);

    return ret;
}

Solution vns(Solution initial, string filePrefix) {
    FILE *objImpFile = nullptr;
    if (!filePrefix.empty())
        objImpFile = fopen(filePrefix.c_str(), "a");

    Solution incumbent = initial;
    Solution delta = convertTo20MhzSol(initial);
    Solution local_min = delta;

    int K_MUL = max(1, n_connections / 100);
    int K_MAX = min(n_connections, 10);
    bool first = true;

    while (!stop()) {
        int k = 1;
        while (!stop() && k <= K_MAX) {
            assert(delta.slots[0].spectrums[3].channels[0].connections.empty());
            delta = local_min;

            perturbation(delta, k * K_MUL);
            computeViolation(delta);

            Solution multiple = multipleRepresentation(delta);

            setDP(multiple);
            delta.violation = calcDP(multiple);

            Solution explicit_sol = local_search(multiple, delta);

            if (essentiallyEqual(explicit_sol.violation, 0.0)) {
                printf("found solution with violation 0.0!\n");
                count_conn(explicit_sol);
                return explicit_sol;
            }

            computeViolation(delta);
            delta = convertTo20MhzSol(explicit_sol);

            if (compareObjectives(delta, local_min) < 0 || first) {
                first = false;
                k = 1;
                local_min = delta;
                rmvUnusedTS(local_min);
            } else
                k += 1;

            if (compareObjectives(local_min, incumbent) < 0 || first) {
                count_conn(explicit_sol);
                double currTime = (((double)(clock() - startTime)) / CLOCKS_PER_SEC);

                if (objImpFile != nullptr)
                    fprintf(objImpFile, "%lf %.2lf\n", currTime, local_min.violation);

                first = false;
                printf("%lf %.2lf %.2lf\n", currTime, local_min.violation, incumbent.violation);
                incumbent = explicit_sol;
                rmvUnusedTS(incumbent);
            }
        }
    }

    if (essentiallyEqual(incumbent.violation, 0.0))
        assert(is_feasible(incumbent));

    return incumbent;
}

Solution reductionHeuristic(char **argv) {
    string objImp = string(argv[3]);
    objImp += "/objective_improvements" + string(argv[2]);
    objImp += ".txt";

    FILE *objImpOut = fopen(objImp.c_str(), "w");
    assert(objImpOut != nullptr);

    Solution S_star = constructive_heuristic(objImp);

    // fprintf(objImpOut, "%.2lf %lu\n", 0.0, S_star.slots.size());
    if (S_star.slots.size() == 1) {
        assert(yFeasible(S_star));
        printf("Constructive heuristic found feasible solution with only one time-slots.\n");
        return S_star;
    }

    printf("initial solution has %lu time-slots\n", S_star.slots.size());
    while (!stop() && S_star.slots.size() > 1) {
        Solution S1 = delete_time_slot(S_star);

        int cnt = 1;
        while (yFeasible(S1) && ++cnt && S1.slots.size() > 1LU) {
            double currTime = (((double)(clock() - startTime)) / CLOCKS_PER_SEC);
            fprintf(objImpOut, "%.3lf %lu\n", currTime, S1.slots.size());

            S_star = S1;
            S1 = delete_time_slot(S1);
        }

        rmvUnusedTS(S1);
        if (yFeasible(S1) && S1.slots.size() == 1LU)
            return S1;

        printf("removed %d ts. Violation is now %lf\n", cnt, S1.violation);
        S1 = vns(S1, string());

        if (yFeasible(S1)) {
            count_conn(S1);
            assert(is_feasible(S1));
            double currTime = (((double)(clock() - startTime)) / CLOCKS_PER_SEC);
            fprintf(objImpOut, "%.3lf %lu\n", currTime, S1.slots.size());

            printf("found sol w/ violation %lf and %lu ts\n", S1.violation, S1.slots.size());
            S_star = S1;
        } else
            puts("vns did not found feasible solution");
    }

    fclose(objImpOut);
    count_conn(S_star);
    return S_star;
}

int main(int argc, char **argv) {
    if (argc != 3 && argc != 5) {
        puts("argument error");
        return 1;
    }

    string path_input = "../instances/md-vrbsp/U_";
    path_input += string(argv[1]) + "/U_";
    // path_input += "/MD-VRBSP_U_";
    path_input += string(argv[1]);
    path_input += "_";
    path_input += string(argv[2]);
    path_input += ".txt";
    freopen(path_input.c_str(), "r", stdin);
    read_data();

    // program_name, instance, version
    if (argc == 3) {
        puts("executing only constructive heuristic and printing results...");
        Solution ch = constructive_heuristic(string()); // TODO
        FILE *aux = fopen("../gurobi/warm.txt", "w");
        print_solution(ch);
        print_solution_to_gurobi(ch, aux);
        return 0;
    }

    maximumTime = stoi(argv[4]) * 1.0;
    Solution inc = reductionHeuristic(argv);
    // print_solution(inc);

    string path_out = string(argv[3]);
    path_out += "/solution" + string(argv[2]);
    path_out += ".txt";
    cout << path_out << endl;
    FILE *file_out = fopen(path_out.c_str(), "w");
    print_solution_to_file(inc, file_out);

    print_solution_to_gurobi(inc);
    return 0;
}
