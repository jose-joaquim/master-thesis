#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <thread>
#include <vector>

#include "gurobi_c++.h"
#include "header.h"

using namespace std;
using mt3 = map<tuple<int, int, int>, GRBVar>;

bool canTransmitUsing(int i, int c) {
  return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]);
}

// -------------- variables --------------------
inline void var_x(GRBModel &model, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsing(i, c)) {
        for (int t = 0; t < T; ++t) {
          string name = "x[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
          x[{i, c, t}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        }
      }
    }
  }
}

inline void var_t(GRBModel &model, vector<GRBVar> &t) {
  for (int i = 0; i < T; ++i) {
    string name = "t[" + to_string(i) + "]";
    t.emplace_back(model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name));
  }
}

inline void var_I(GRBModel &model, vector<GRBVar> &I) {
  for (int i = 0; i < N; ++i) {
    string name = "I[" + to_string(i) + "]";
    I.emplace_back(model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name));
  }
}

inline void var_Iij(GRBModel &model, mt3 &Iij) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsing(i, c)) {
        for (int t = 0; t < T; ++t) {
          string name = "Iij[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
          Iij[{i, c, t}] =
              model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
        }
      }
    }
  }
}

inline void var_z(GRBModel &model, mt3 &z) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsing(i, c)) {
        for (int t = 0; t < T; ++t) {
          string name = "z[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
          z[{i, c, t}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        }
      }
    }
  }
}

// -------------- constraints ------------------

inline void symmetry1(GRBModel &model, vector<GRBVar> &t) {
  for (int i = 0; i < T - 1; ++i) {
    string name = "ts" + to_string(i);
    model.addConstr(t[i + 1] <= t[i], name);
  }
}

inline void symmetry2(GRBModel &model, mt3 &x) {
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
        model.addConstr(expr == 0, name);
      }
    }
  }
}

inline void unique(GRBModel &model, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c))
        continue;

      for (int t = 0; t < T; ++t)
        expr += x[{i, c, t}];
    }

    string name = "unique[" + to_string(i) + "]";
    model.addConstr(expr == 1, name);
  }
}

inline void waste(GRBModel &model, mt3 &x, vector<GRBVar> &t) {
  for (int i = 0; i < N; ++i) {
    for (int _t = 0; _t < T; ++_t) {
      GRBLinExpr expr = 0;
      for (int c = 0; c < C; ++c)
        if (canTransmitUsing(i, c))
          expr += x[{i, c, _t}];

      string name = "waste[" + to_string(i) + "," + to_string(_t) + "]";
      model.addConstr(expr <= t[_t], name);
    }
  }
}

inline void ch_overlap(GRBModel &model, mt3 &z, mt3 &x) {
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

        string name = "over[" + to_string(i) + "," + to_string(c1) + "," +
                      to_string(t) + "]";
        model.addConstr(expr == z[{i, c1, t}], name);
      }
    }
  }
}

inline void interch(GRBModel &model, mt3 &z, mt3 &Iij) {
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

        string name = "interch[" + to_string(i) + "," + to_string(c) + "," +
                      to_string(t) + "]";
        model.addConstr(expr == Iij[{i, c, t}], name);
      }
    }
  }
}

inline void bigG(GRBModel &model, vector<GRBVar> &I, mt3 &Iij, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c))
        continue;

      for (int t = 0; t < T; ++t) {
        string name = "bigG[" + to_string(i) + "," + to_string(c) + "," +
                      to_string(t) + "]";
        model.addConstr(I[i] >= Iij[{i, c, t}] - BM[i] * (1 - x[{i, c, t}]),
                        name);
      }
    }
  }
}

inline void bigL(GRBModel &model, vector<GRBVar> &I, mt3 &Iij, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c))
        continue;

      for (int t = 0; t < T; ++t) {
        string name = "bigL[" + to_string(i) + "," + to_string(c) + "," +
                      to_string(t) + "]";
        model.addConstr(I[i] <= Iij[{i, c, t}] + BM[i] * (1 - x[{i, c, t}]),
                        name);
      }
    }
  }
}

inline void sinr(GRBModel &model, vector<GRBVar> &I, mt3 &x) {
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
    model.addConstr(I[i] <= expr, name);
  }
}

void build_model(GRBModel &model, mt3 &x, mt3 &z, vector<GRBVar> &I, mt3 &Iij,
                 vector<GRBVar> &t) {
  // variables
  // printf("variables...\n");
  var_x(model, x);
  var_z(model, z);
  var_Iij(model, Iij);
  var_I(model, I);
  var_t(model, t);

  // constraints
  model.update();
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
}

int main(int argc, char *argv[]) {
  read_data();

  try {
    GRBEnv env;
    GRBModel master(env);

    // L_0: |L| time-slots, cada um com os links alocados no canal de 160 MHz
    vector<GRBVar> vars;
    for (int l = 0; l < N; ++l)
      vars.emplace_back(
          master.addVar(0.0, 1.0, 1.0, GRB_BINARY, "sol0_" + to_string(l)));

    // sum_{l \in L} \sum_{c \in C} \sum_{t \in T} x_{l}^{ct} >= 1
    vector<GRBConstr> masterConstrs;
    for (int l = 0; l < N; ++l)
      masterConstrs.emplace_back(master.addConstr(vars[l] >= 1.0));

    // master.optimize();
    int loops = 0;
    while (true) {
      loops++;
      master.update();
      master.write("master" + to_string(loops) + ".lp");

      GRBModel relax = master.relax();
      relax.set(GRB_DoubleParam_IntFeasTol, 1e-5);
      relax.set(GRB_IntParam_LogToConsole, 0);
      relax.optimize();

      cout << "relaxed master solved! objective: "
           << relax.get(GRB_DoubleAttr_ObjVal) << " dual variables:" << endl;
      vector<double> pi;
      for (int l = 0; l < N; ++l) {
        pi.emplace_back(relax.getConstr(l).get(GRB_DoubleAttr_Pi));
        cout << pi.back() << " ";
      }

      cout << endl;
      GRBModel sub(env);
      vector<GRBVar> I, t;
      mt3 x, Iij, z;
      build_model(sub, x, z, I, Iij, t);

      GRBLinExpr expr = 0;
      for (int l = 0; l < N; ++l)
        for (int c = 0; c < C; ++c)
          if (canTransmitUsing(l, c)) {
            for (int t = 0; t < T; ++t)
              expr += x[{l, c, t}] * pi[l];
          }

      sub.setObjective(expr, GRB_MAXIMIZE);

      sub.set(GRB_DoubleParam_IntFeasTol, 1e-5);
      sub.update();
      sub.write("sub" + to_string(loops) + ".lp");
      sub.set(GRB_IntParam_LogToConsole, 0);
      sub.optimize();

      if (sub.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        cout << "sub solved! objective: " << sub.get(GRB_DoubleAttr_ObjVal)
             << endl;

      map<tuple<int, int, int>, int> solution;
      for (int l = 0; l < N; ++l)
        for (int c = 0; c < C; ++c)
          if (canTransmitUsing(l, c)) {
            for (int t = 0; t < T; ++t) {
              solution[{l, c, t}] =
                  static_cast<int>(x[{l, c, t}].get(GRB_DoubleAttr_X) + .5);
            }
          }

      map<int, GRBColumn> col;
      for (int l = 0; l < N; ++l)
        for (int c = 0; c < C; ++c)
          if (canTransmitUsing(l, c)) {
            for (int t = 0; t < T; ++t)
              if (solution[{l, c, t}] > 0) {
                // cout << l << " " << c << " " << t << endl;
                col[t].addTerm(solution[{l, c, t}], masterConstrs[l]);
              }
          }

      for (auto &e : col) {
        vars.emplace_back(
            master.addVar(0.0, 1.0, 1.0, GRB_BINARY, e.second,
                          "x" + to_string(loops) + "_" + to_string(e.first)));
      }

      double stop = 1 - sub.get(GRB_DoubleAttr_ObjVal);
      if (approximatelyEqual(stop, 0.0) || definitelyLessThan(stop, 0.0))
        break;
    }

    cout << "======================" << endl;
    cout << "solving last model" << endl;
    master.update();
    master.write("masterFinal.lp");
    master.set(GRB_IntParam_LogToConsole, 0);
    master.set(GRB_DoubleParam_IntFeasTol, 1e-5);
    master.optimize();

    FILE* obj = fopen("obj", "a");
    fprintf(obj, "%lf\n", master.get(GRB_DoubleAttr_ObjVal));
    fclose(obj);
    cout <<  master.get(GRB_DoubleAttr_ObjVal) << endl;
  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
