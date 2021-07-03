#include "VNS.h"

bool is_feasible(const Solution &ret) {
    int cnt = 0;
    for (const TimeSlot &ts : ret.slots) {
        for (const Spectrum &sp : ts.spectrums) {
            for (const Channel &ch : sp.channels) {
                for (const Connection &conn : ch.connections) {
                    // printf("conn %d\n", conn.id);
                    cnt += 1;
                    if (conn.throughput < gma[conn.id]) {
                        printf("%lf %lf\n", conn.throughput, gma[conn.id]);
                        return false;
                    }
                }
            }
        }
    }

    // printf("%d %d\n", cnt, n_connections);
    return cnt == n_connections;
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

    for (Connection &conn : ch.connections) {
        if (conn.throughput < gma[conn.id])
            return false;
    }

    return true;
}

Solution constructive_heuristic() {
    TimeSlot dummy_ts(init_conf);
    Solution ret({dummy_ts});

    vector<int> links;
    for (int i = 0; i < n_connections; i++)
        links.emplace_back(i);

    shuffle(links.begin(), links.end(), whatever);
    for (int conn : links) {
        // printf("tentando inserir conexao %d\n", conn);
        bool success = false;
        for (int t = 0; t < ret.slots.size(); t++) {
            for (int s = 0; s < ret.slots[t].spectrums.size(); s++) {
                if (success)
                    break;

                for (int c = 0; c < ret.slots[t].spectrums[s].channels.size(); c++) {
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

    computeThroughput(ret);
    assert(is_feasible(ret));
    return ret;
}

void delete_time_slot(Solution &sol) {
    int to_del = rng.randInt(sol.slots.size() - 1);
    vector<int> links_id;
    for (const Spectrum &sp : sol.slots[to_del].spectrums) {
        for (const Channel &ch : sp.channels) {
            for (const Connection &conn : ch.connections) {
                links_id.push_back(conn.id);
            }
        }
    }

    swap(sol.slots[to_del], sol.slots.back());
    sol.slots.pop_back();

    sol.unscheduled = links_id;
}

Solution vns(string filePrefix) {
    // 1st step: constructive heuristic
    // 2nd step: dp algorithm
    // 3rd step: local search
    // 4th step: main loop until stopping condition
    //         : 1) perturbation
    //         : 2) local search
    //         : 3) update

    FILE *file_comparative = nullptr;
    if (!filePrefix.empty())
        file_comparative = fopen(filePrefix.c_str(), "w");

    Solution init_sol = constructive_heuristic(); // TODO (?)
    compute_violation(init_sol);

    // if (essentiallyEqual(init_sol.violation, 0.0)) {
    //     puts("opa");
    //     return init_sol;
    // }

    Solution delta = convertTo20MhzSol(init_sol); // DONE
    Solution rep = multipleRepresentation(delta); // DONE
    setDP(rep);
    double retOF = calcDP(rep);

    // Then, reconstruct optimal local solution
    Solution incumbent = reconstruct_sol(rep); // DONE
    incumbent.throughput = init_sol.throughput;
    fprintf(file_comparative, "%lf %u %lf\n", 0.0, incumbent.slots.size(), incumbent.violation);

    Solution local_min = delta;

    int K_MUL = max(1, n_connections / 100);
    int K_MAX = 10;
    bool first = true;
    startTime = clock();
    while (!stop()) {
        int k = 1;
        while (!stop() && k <= K_MAX) {
            delta = local_min;

            if (delta.slots.size() > 1)
                delete_time_slot(delta);

            perturbation(delta, k * K_MUL); // DONE
            compute_violation(delta);

            Solution multiple = multipleRepresentation(delta); // DONE

            setDP(multiple);
            delta.throughput = calcDP(multiple);

            // printf("%lf %lf\n", incumbent.throughput, delta.throughput);

            Solution explicit_sol = local_search(multiple, delta); // TODO
            fix_channels(explicit_sol);                            // DONE

            compute_violation(delta);
            delta = convertTo20MhzSol(explicit_sol); // DONE

            if (essentiallyEqual(delta.violation, 0.0)) {
                puts("hm 2");
                return delta;
            }

            // cout << delta.violation << " " << local_min.violation << endl;
            if (definitelyLessThan(delta.violation, local_min.violation) || first) {
                first = false;
                k = 1;
                local_min = delta;
            } else {
                k += 1;
            }
            if (definitelyLessThan(local_min.violation, incumbent.violation) || first) {
                double elapsed_time = (((double)(clock() - startTime)) / CLOCKS_PER_SEC);
                if (file_comparative != nullptr)
                    fprintf(file_comparative, "%lf %u %lf\n", elapsed_time, incumbent.slots.size(), local_min.violation);
                
                first = false;
                printf("melhorei %lf %.3lf %.3lf\n", elapsed_time, local_min.violation, incumbent.violation);
                incumbent = explicit_sol;
            }
        }
    }

    return incumbent;
}

int main(int argc, char **argv) {
    if (argc != 3 && argc != 5) {
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
    read_data();

    // program_name, instance, version
    if (argc == 3) {
        puts("executing only constructive heuristic and printing results...");
        Solution ch = constructive_heuristic();
        FILE *aux = fopen("../gurobi/warm.txt", "w");
        print_solution(ch);
        print_solution_to_gurobi(ch, aux);
        return 0;
    }

    maximumTime = stoi(argv[4]) * 1.0;
    Solution inc = vns();
    print_solution(inc);
    return 0;
}
