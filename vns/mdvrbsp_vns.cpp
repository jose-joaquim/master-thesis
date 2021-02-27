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

Solution vns() {
    // 1st step: constructive heuristic
    // 2nd step: dp algorithm
    // 3rd step: local search
    // 4th step: main loop until stopping condition
    //         : 1) perturbation
    //         : 2) local search
    //         : 3) update

    Solution init_sol = constructive_heuristic(); // TODO (?)
    Solution delta = convertTo20MhzSol(init_sol); // DONE

    // // Find best partitioning
    // Solution multi_rep = multipleRepresentation(delta); // DONE
    // runDP(multi_rep);                                      // DONE
    //
    // // Then, reconstruct optimal local solution
    // Solution global_max = reconstruct_sol(multi_rep); // DONE
    // Solution local_max = delta;
    //
    // int K_MUL = max(1, n_connections / 100);
    // int K_MAX = 10;
    // startTime = clock();
    // while (!stop()) {
    //     int k = 1;
    //     while (!stop() && k <= K_MAX) {
    //         delta = local_min;
    //
    //         delete_time_slot(delta);
    //         perturbation(delta, k * K_MUL);            // DONE
    //         multi_rep = multipleRepresentation(delta); // DONE
    //         runDP(multi_rep);                             // DONE
    //
    //         Solution explicit_sol = local_search(multi_rep, delta); // TODO
    //         fix_channels(explicit_sol);                                // DONE
    //
    //         delta = convertTo20MhzSol(explicit_sol); // DONE
    //
    //         if (delta < local_min) {
    //             k = 1;
    //             local_min = delta;
    //         } else {
    //             k += 1;
    //         }
    //         if (local_min < global_min) {
    //             global_min = explicit_sol;
    //         }
    //     }
    // }
    //
    // return global_min;
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
    read_data(); // Input will be redirected in the command call
    vns();
    return 0;
}
