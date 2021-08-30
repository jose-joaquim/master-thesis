// run 16 1 results/vns-vrbsp/U_16 10
// URGENTE: verificar se quando o canal eh inviavel se todas as conexoes tem throughput setado para
// 0!!!!
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

    if (retOption)
        return cnt == n_connections;

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
                    if (t == 0 && s == 3 && c == 0)
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
    assert(is_feasible(ret));
    return ret;
}

Solution delete_time_slot(const Solution &sol) {
    Solution ret = Solution(sol);

    int to_del = rng.randInt(ret.slots.size() - 1);
    vector<Connection> links;
    for (const Spectrum &sp : ret.slots[to_del].spectrums) {
        for (const Channel &ch : sp.channels) {
            for (const Connection &conn : ch.connections) {
                Connection aux(conn);
                aux.throughput = 0.0;
                aux.interference = 0.0;
                aux.SINR = 0.0;
                links.emplace_back(aux);
            }
        }
    }

    swap(ret.slots[to_del], ret.slots.back());
    ret.slots.pop_back();
    ret.slots[0].spectrums[3].channels[0].connections = links;
    K_AddDrop(ret, links.size());
    computeViolation(ret);

    count_conn(ret);
    return ret;
}

Solution vns(Solution initial, string filePrefix) {
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

    count_conn(initial);
    Solution delta = convertTo20MhzSol(initial);
    count_conn(delta);
    Solution rep = multipleRepresentation(delta);
    setDP(rep);
    double retOF = calcDP(rep);

    // Then, reconstruct optimal local solution
    Solution incumbent = reconstruct_sol(rep);
    incumbent.throughput = initial.throughput;

    if (file_comparative != nullptr)
        fprintf(file_comparative, "%lf %lu %lf\n", 0.0, incumbent.slots.size(),
                incumbent.violation);

    Solution local_min = delta;
    count_conn(local_min);
    
    int K_MUL = max(1, n_connections / 100);
    int K_MAX = 10;
    bool first = true;
    startTime = clock();
    while (!stop()) {
        int k = 1;
        while (!stop() && k <= K_MAX) {
            puts("begin");
            
            delta = local_min;

            count_conn(delta);
            perturbation(delta, k * K_MUL);

            count_conn(delta);
            computeViolation(delta);

            Solution multiple = multipleRepresentation(delta);
            
            count_conn(multiple, true);
            setDP(multiple);
            delta.violation = calcDP(multiple);
            
            Solution explicit_sol = local_search(multiple, delta);
            count_conn(explicit_sol);
            
            if (essentiallyEqual(explicit_sol.violation, 0.0)) {
                printf("found solution with violation 0.0!\n");
                return explicit_sol;
            }
            
            fix_channels(explicit_sol);
            
            computeViolation(delta);
            delta = convertTo20MhzSol(explicit_sol);
            count_conn(delta);
            
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
                    fprintf(file_comparative, "%lf %lu %lf\n", elapsed_time, incumbent.slots.size(), local_min.violation);
            
                // assert(is_feasible(explicit_sol));
                first = false;
                printf("melhorei %lf %.3lf %.3lf\n", elapsed_time, local_min.violation,
                       incumbent.violation);
                incumbent = explicit_sol;
            }
            
            count_conn(local_min);
            count_conn(delta);
            count_conn(explicit_sol);
            count_conn(incumbent);
            // puts("end iteration");
        }
    }

    count_conn(incumbent);

    return incumbent;
}

Solution reductionHeuristic() {
    Solution S_star = constructive_heuristic();
    printf("initial solution has %lu time-slots\n", S_star.slots.size());
    while (not stop()) {
        Solution S1 = delete_time_slot(S_star);

        int cnt = 1;
        while (yFeasible(S1) && ++cnt)
            S1 = delete_time_slot(S1);
        
        printf("removed %d time-slots. Violation is now %lf\n", cnt, S1.violation);

        S1 = vns(S1, string());

        puts("leaving vns...");
        if (yFeasible(S1)) {
            printf("found feasible solution with violation %lf and %lu time-slots\n", S1.violation,
                   S1.slots.size());
            assert(is_feasible(S1));
            S_star = S1;
        }
    }

    return S_star;
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
    Solution inc = reductionHeuristic();
    print_solution(inc);
    return 0;
}
