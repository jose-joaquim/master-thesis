#include "VNS.h"

using namespace std;

bool is_feasible(FILE *file) {
    // TODO: no arquivo de solucao preciso dizer a largura de banda para o computo correto do data-rate
    Solution sol;
    double global_throughput;
    int n_slots;
    fscanf(file, "%lf %d", &global_throughput, &n_slots);
    printf("li %lf %d\n", global_throughput, n_slots);
    vector<TimeSlot> v_slots;
    for (int i = 0; i < n_slots; i++) {
        TimeSlot new_ts;
        int n_specs;
        fscanf(file, "%d", &n_specs);
        printf("= 1 = li %d specs\n", n_specs);
        vector<Spectrum> v_specs(n_specs);
        for (int j = 0; j < n_specs; j++) {
            int n_channels;
            fscanf(file, "%d", &n_channels);
            printf("= 2 = li %d canais\n", n_channels);
            Spectrum new_spec;
            vector<Channel> v_channels(n_channels);
            for (int k = 0; k < n_channels; k++) {
                int n_conns, bandwidth;
                fscanf(file, "%d %d", &bandwidth, &n_conns);
                Channel new_channel(bandwidth);
                printf("= 3 = li %d conns\n", n_conns);
                vector<Connection> conns;
                for (int l = 0; l < n_conns; l++) {
                    int indice;
                    fscanf(file, "%d", &indice);
                    printf("inseri conn %d\n", indice);
                    conns.emplace_back(indice);
                }

                new_channel.connections = conns;
                v_channels.emplace_back(new_channel);
            }
            new_spec.channels = v_channels;
            v_specs.emplace_back(new_spec);
        }

        new_ts.spectrums = v_specs;
        v_slots.emplace_back(new_ts);
    }

    sol.slots = v_slots;
    computeThroughput(sol);

    printf("%lf\n", sol.throughput);
    
#ifdef MDVRBSP
    int cnt = 0;
    for (const TimeSlot &ts : sol.slots) {
        for (const Spectrum &sp : ts.spectrums) {
            for (const Channel &ch : sp.channels) {
                for (const Connection &conn : ch.connections) {
                    cnt += 1;
                    if (conn.throughput < gma[conn.id]) {
                        printf("%lf %lf\n", conn.throughput, gma[conn.id]);
                        return false;
                    }
                }
            }
        }
    }

    if (cnt != n_connections) {
        // printf("(2) %d %d\n", cnt, n_connections);
        return false;
    }
#else
    if (!double_equals(sol.throughput, global_throughput)) {
        // puts("opaaaaaa");
        return false;
    }
#endif
    return true;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        puts("argument error in feasibility check");
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

    FILE *sol_file = fopen(argv[3], "r");
    bool ans = is_feasible(sol_file);
    printf("%s\n", ans ? "TUDO CERTO" : "DEU ERRADO");

    fclose(sol_file);
    return 0;
}
