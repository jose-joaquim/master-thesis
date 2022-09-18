#include "header.h"
#include "gurobi_c++.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>

using namespace std;

// -------------- variables --------------------
inline void var_x(GRBModel *model, GRBVar ***&x) {
    x = new GRBVar **[N];
    for (int i = 0; i < N; ++i) {
        x[i] = new GRBVar *[C];
        for (int c = 0; c < C; ++c) {
            x[i][c] = new GRBVar[T];
            for (int t = 0; t < T; ++t) {
                // if (definitelyGreaterThan(GMM[i], DR[11][cToBIdx(c)])) {
                string name = "x[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                x[i][c][t] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                //}
            }
        }
    }
}

inline void var_t(GRBModel *model, GRBVar *&t) {
    t = new GRBVar[T];
    for (int i = 0; i < T; ++i) {
        string name = "t[" + to_string(i) + "]";
        t[i] = model->addVar(0.0, 1.0, 1.0, GRB_BINARY, name);
    }
}

inline void var_I(GRBModel *model, GRBVar *&I) {
    I = new GRBVar[N];
    for (int i = 0; i < N; ++i) {
        string name = "I[" + to_string(i) + "]";
        I[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
    }
}

inline void var_Iij(GRBModel *model, GRBVar ***&Iij) {
    Iij = new GRBVar **[N];
    for (int i = 0; i < N; ++i) {
        Iij[i] = new GRBVar *[C];
        for (int c = 0; c < C; ++c) {
            Iij[i][c] = new GRBVar[T];
            for (int t = 0; t < T; ++t) {
                string name = "Iij[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                Iij[i][c][t] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
            }
        }
    }
}

inline void var_z(GRBModel *model, GRBVar ***&z) {
    z = new GRBVar **[N];
    for (int i = 0; i < N; ++i) {
        z[i] = new GRBVar *[C];
        for (int c = 0; c < C; ++c) {
            z[i][c] = new GRBVar[T];
            for (int t = 0; t < T; ++t) {
                string name = "z[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                z[i][c][t] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
            }
        }
    }
}

// -------------- constraints ------------------

inline void symmetry1(GRBModel *model, GRBVar *t) {
    for (int i = 0; i < T - 1; ++i) {
        string name = "ts" + to_string(i);
        model->addConstr(t[i + 1] <= t[i], name);
    }
}

inline void symmetry2(GRBModel *model, GRBVar ***x) {
    for (int t = 0; t < T; ++t) {
        for (int i = 0; i < N; ++i) {
            if (i < t) {
                GRBLinExpr expr = 0;
                for (int c = 0; c < C; ++c) {
                    expr += x[i][c][t];
                }

                string name = "swap[" + to_string(t) + "," + to_string(i) + "]";
                model->addConstr(expr == 0, name);
            }
        }
    }
}

inline void unique(GRBModel *model, GRBVar ***&x) {
    for (int i = 0; i < N; ++i) {
        GRBLinExpr expr = 0;
        for (int c = 0; c < C; ++c) {
            for (int t = 0; t < T; ++t) {
                expr += x[i][c][t];
            }
        }

        string name = "unique" + to_string(i);
        model->addConstr(expr == 1, name);
    }
}

inline void waste(GRBModel *model, GRBVar ***x, GRBVar *t) {
    for (int i = 0; i < N; ++i) {
        for (int _t = 0; _t < T; ++_t) {
            GRBLinExpr expr = 0;
            for (int c = 0; c < C; ++c) {
                expr += x[i][c][_t];
            }

            string name = "waste[" + to_string(i) + "," + to_string(_t) + "]";
            model->addConstr(expr <= t[_t], name);
        }
    }
}

inline void ch_overlap(GRBModel *model, GRBVar ***z, GRBVar ***x) {
    for (int i = 0; i < N; ++i) {
        for (int t = 0; t < T; ++t) {
            for (int c1 = 0; c1 < C; ++c1) {
                GRBLinExpr expr = 0;
                for (int c2 = 0; c2 < C; ++c2) {
                    if (overlap[c1][c2])
                        expr += x[i][c2][t];
                }

                string name =
                    "over[" + to_string(i) + "," + to_string(t) + "," + to_string(c1) + "]";
                model->addConstr(expr == z[i][c1][t], name);
            }
        }
    }
}

inline void interch(GRBModel *model, GRBVar ***z, GRBVar ***Iij) {
    for (int i = 0; i < N; ++i) {
        for (int t = 0; t < T; ++t) {
            for (int c = 0; c < C; ++c) {
                GRBLinExpr expr = 0;
                for (int u = 0; u < N; ++u) {
                    if (u != i)
                        expr += AFF[u][i] * z[u][c][t];
                }

                string name =
                    "interch[" + to_string(i) + "," + to_string(t) + "," + to_string(c) + "]";
                model->addConstr(expr == Iij[i][c][t], name);
            }
        }
    }
}

inline void bigG(GRBModel *model, GRBVar *I, GRBVar ***Iij, GRBVar ***x) {
    for (int i = 0; i < N; ++i) {
        for (int t = 0; t < T; ++t) {
            for (int c = 0; c < C; ++c) {
                string name =
                    "bigG[" + to_string(i) + "," + to_string(t) + "," + to_string(c) + "]";
                model->addConstr(I[i] >= Iij[i][c][t] - BM[i] * (1 - x[i][c][t]), name);
            }
        }
    }
}

inline void bigL(GRBModel *model, GRBVar *I, GRBVar ***Iij, GRBVar ***x) {
    for (int i = 0; i < N; ++i) {
        for (int t = 0; t < T; ++t) {
            for (int c = 0; c < C; ++c) {
                string name =
                    "bigL[" + to_string(i) + "," + to_string(t) + "," + to_string(c) + "]";
                model->addConstr(I[i] <= Iij[i][c][t] + BM[i] * (1 - x[i][c][t]), name);
            }
        }
    }
}

inline void sinr(GRBModel *model, GRBVar *I, GRBVar ***x) {
    for (int i = 0; i < N; ++i) {
        GRBLinExpr expr = 0;
        for (int t = 0; t < T; ++t) {
            for (int c = 0; c < C; ++c) {
                expr += (AFF[i][i] / B[i][cToBIdx(c)] - NOI) * x[i][c][t];
            }
        }

        string name = "sinr[" + to_string(i) + "]";
        model->addConstr(I[i] <= expr, name);
    }
}

inline void cant(GRBModel *model, GRBVar ***x) {
    for (int c = 0; c < C; ++c) {
        for (int i = 0; i < N; ++i) {
            for (int t = 0; t < T; ++t) {
                if (DR[11][cToBIdx(c)] < GMM[i]) {
                    string name =
                        "cant[" + to_string(c) + "," + to_string(i) + "," + to_string(t) + "]";
                    model->addConstr(x[i][c][t] == 0, name);
                }
            }
        }
    }
}

// ---------------------------------------------

int main(int argc, char **argv) {
    read_data();
    printf("%d %d %d\n", N, C, T);
    auto start = high_resolution_clock::now();
    GRBEnv env;
    GRBModel *model = new GRBModel(env);
    GRBVar ***x, ***Iij, ***z, *I, *t;
    // variables
    // printf("variables...\n");
    var_x(model, x);
    var_I(model, I);
    var_z(model, z);
    var_t(model, t);
    var_Iij(model, Iij);

    // constraints
    model->update();
    // printf("constraints...\n");
    // cant(model, x);
    symmetry1(model, t);
    symmetry2(model, x);
    unique(model, x);
    waste(model, x, t);
    ch_overlap(model, z, x);
    interch(model, z, Iij);
    bigG(model, I, Iij, x);
    bigL(model, I, Iij, x);
    sinr(model, I, x);
    auto stop = high_resolution_clock::now();
    duration<double> ms_double = stop - start;
    cout << ms_double.count() << endl;
    // optimize
    model->update();
    model->write("seila.lp");
    model->set(GRB_IntParam_LogToConsole, 0);
    model->optimize();

    if (model->get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        FILE *result = fopen("result", "a");
        fprintf(result, "%lf\n", model->get(GRB_DoubleAttr_ObjVal));
        fclose(result);
    } else
        cout << "not optimal!" << endl;
    return 0;
}
