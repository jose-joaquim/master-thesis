#include "header.h"
#include "gurobi_c++.h"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <tuple>

using namespace std;
using mt3 = map<tuple<int, int, int>, GRBVar>;

bool canTransmitUsing(int i, int c) { return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]); }

// -------------- variables --------------------
inline void var_x(GRBModel *model, mt3 &x) {
    for (int i = 0; i < N; ++i) {
        for (int c = 0; c < C; ++c) {
            if (canTransmitUsing(i, c)) {
                for (int t = 0; t < T; ++t) {
                    string name =
                        "x[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                    x[{i, c, t}] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                }
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

inline void var_Iij(GRBModel *model, mt3 &Iij) {
    for (int i = 0; i < N; ++i) {
        for (int c = 0; c < C; ++c) {
            if (canTransmitUsing(i, c)) {
                for (int t = 0; t < T; ++t) {
                    string name =
                        "Iij[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                    Iij[{i, c, t}] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
                }
            }
        }
    }
}

inline void var_z(GRBModel *model, mt3 &z) {
    for (int i = 0; i < N; ++i) {
        for (int c = 0; c < C; ++c) {
            if (canTransmitUsing(i, c)) {
                for (int t = 0; t < T; ++t) {
                    string name =
                        "z[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                    z[{i, c, t}] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                }
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

inline void symmetry2(GRBModel *model, mt3 &x) {
    for (int t = 0; t < T; ++t) {
        for (int i = 0; i < N; ++i) {
            if (i < t) {
                GRBLinExpr expr = 0;
                for (int c = 0; c < C; ++c) {
                    if (canTransmitUsing(i, c)) {
                        expr += x[{i, c, t}];
                    }
                }

                string name = "swap[" + to_string(t) + "," + to_string(i) + "]";
                model->addConstr(expr == 0, name);
            }
        }
    }
}

inline void unique(GRBModel *model, mt3 &x) {
    for (int i = 0; i < N; ++i) {
        GRBLinExpr expr = 0;
        for (int c = 0; c < C; ++c) {
            if (!canTransmitUsing(i, c))
                continue;

            for (int t = 0; t < T; ++t)
                expr += x[{i, c, t}];
        }

        string name = "unique[" + to_string(i) + "]";
        model->addConstr(expr == 1, name);
    }
}

inline void waste(GRBModel *model, mt3 &x, GRBVar *t) {
    for (int i = 0; i < N; ++i) {
        for (int _t = 0; _t < T; ++_t) {
            GRBLinExpr expr = 0;
            for (int c = 0; c < C; ++c)
                if (canTransmitUsing(i, c))
                    expr += x[{i, c, _t}];

            string name = "waste[" + to_string(i) + "," + to_string(_t) + "]";
            model->addConstr(expr <= t[_t], name);
        }
    }
}

inline void ch_overlap(GRBModel *model, mt3 &z, mt3 &x) {
    for (int i = 0; i < N; ++i) {
        for (int c1 = 0; c1 < C; ++c1) {
            if (!canTransmitUsing(i, c1))
                continue;

            for (int t = 0; t < T; ++t) {
                GRBLinExpr expr = 0;
                for (int c2 = 0; c2 < C; ++c2) {
                    if (overlap[c1][c2] && canTransmitUsing(i, c2))
                        expr += x[{i, c2, t}];
                }

                string name =
                    "over[" + to_string(i) + "," + to_string(c1) + "," + to_string(t) + "]";
                model->addConstr(expr == z[{i, c1, t}], name);
            }
        }
    }
}

inline void interch(GRBModel *model, mt3 &z, mt3 &Iij) {
    for (int i = 0; i < N; ++i) {
        for (int t = 0; t < T; ++t) {
            for (int c = 0; c < C; ++c) {
                if (!canTransmitUsing(i, c))
                    continue;

                GRBLinExpr expr = 0;
                for (int u = 0; u < N; ++u) {
                    if (u != i && canTransmitUsing(u, c))
                        expr += AFF[u][i] * z[{u, c, t}];
                }

                string name =
                    "interch[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                model->addConstr(expr == Iij[{i, c, t}], name);
            }
        }
    }
}

inline void bigG(GRBModel *model, GRBVar *I, mt3 &Iij, mt3 &x) {
    for (int i = 0; i < N; ++i) {
        for (int c = 0; c < C; ++c) {
            if (!canTransmitUsing(i, c))
                continue;

            for (int t = 0; t < T; ++t) {
                string name =
                    "bigG[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                model->addConstr(I[i] >= Iij[{i, c, t}] - BM[i] * (1 - x[{i, c, t}]), name);
            }
        }
    }
}

inline void bigL(GRBModel *model, GRBVar *I, mt3 &Iij, mt3 &x) {
    for (int i = 0; i < N; ++i) {
        for (int c = 0; c < C; ++c) {
            if (!canTransmitUsing(i, c))
                continue;

            for (int t = 0; t < T; ++t) {
                string name =
                    "bigL[" + to_string(i) + "," + to_string(c) + "," + to_string(t) + "]";
                model->addConstr(I[i] <= Iij[{i, c, t}] + BM[i] * (1 - x[{i, c, t}]), name);
            }
        }
    }
}

inline void sinr(GRBModel *model, GRBVar *I, mt3 &x) {
    for (int i = 0; i < N; ++i) {
        GRBLinExpr expr = 0;
        for (int c = 0; c < C; ++c) {
            if (!canTransmitUsing(i, c))
                continue;

            for (int t = 0; t < T; ++t) {
                expr += (AFF[i][i] / B[i][cToBIdx(c)] - NOI) * x[{i, c, t}];
            }
        }

        string name = "sinr[" + to_string(i) + "]";
        model->addConstr(I[i] <= expr, name);
    }
}

// ---------------------------------------------

int main(int argc, char **argv) {
    read_data();
    // auto start = high_resolution_clock::now();
    GRBEnv env;
    GRBModel *model = new GRBModel(env);
    GRBVar *I, *t;
    mt3 x, Iij, z;
    // variables
    // printf("variables...\n");
    var_x(model, x);
    var_z(model, z);
    var_Iij(model, Iij);
    var_I(model, I);
    var_t(model, t);

    // constraints
    model->update();
    // printf("constraints...\n");
    symmetry1(model, t);
    symmetry2(model, x);
    unique(model, x);
    waste(model, x, t);
    ch_overlap(model, z, x);
    interch(model, z, Iij);
    bigG(model, I, Iij, x);
    bigL(model, I, Iij, x);
    sinr(model, I, x);
    // auto stop = high_resolution_clock::now();
    // duration<double> ms_double = stop - start;
    // cout << ms_double.count() << endl;
    // optimize
    model->update();
    // model->write("seila.lp");
    model->set(GRB_IntParam_LogToConsole, 0);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->optimize();
    
    if (model->get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
        // model->write("sol.sol");
        string res = "result" + string(argv[1]);
        FILE *result = fopen(res.c_str(), "a");
        fprintf(result,
                "%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n",
                argv[2],
                model->get(GRB_DoubleAttr_ObjVal),
                model->get(GRB_DoubleAttr_ObjBoundC),
                model->get(GRB_DoubleAttr_MIPGap),
                model->get(GRB_IntAttr_NumVars),
                model->get(GRB_IntAttr_NumConstrs),
                model->get(GRB_DoubleAttr_DNumNZs),
                model->get(GRB_DoubleAttr_NodeCount),
                model->get(GRB_DoubleAttr_Runtime));
        fclose(result);
    } else {
        model->set(GRB_IntParam_IISMethod, 1);
        model->computeIIS();
        model->write("xi.ilp");
        cout << "not optimal!" << endl;
    }
    return 0;
}
