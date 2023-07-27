#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <tuple>
#include <type_traits>

#include "gurobi_c++.h"
#include "header.h"

using namespace std;
using mt3 = map<tuple<int, int, int>, GRBVar>;
using mt2 = map<tuple<int, int>, GRBVar>;

vector<vector<int>> C_b = {
    {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
        13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
    },
    {25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36},
    {37, 38, 39, 40, 41, 42},
    {43, 44},
};

bool canTransmitUsingChannel(int i, int c) {
  return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]) ||
         approximatelyEqual(DR[11][cToBIdx(c)], GMM[i]);
}

bool canTransmitUsingBandwidth(int i, int b, int m) {
  return definitelyGreaterThan(DR[m][b], GMM[i]) ||
         approximatelyEqual(DR[m][b], GMM[i]);
}

// -------------- variables --------------------
inline void var_y(GRBModel *model, mt3 &y, vector<int> &conns) {
  int cnt = 0;
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];

    for (int b = 0; b < 4; ++b) {
      for (int m = 0; m < 12; ++m) {
        if (canTransmitUsingBandwidth(i, b, m)) {
          string name = "y[" + to_string(i) + "," + to_string(b) + "," +
                        to_string(m) + "]";
          y[{i, b, m}] = model->addVar(0.0, 1.0, DR[m][b], GRB_BINARY, name);
          cnt++;
        }
      }
    }
  }
  cout << "created " << cnt << " y variables" << endl;
}

inline void var_x(GRBModel *model, mt2 &x, vector<int> &conns) {
  int cnt = 0;
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsingChannel(i, c)) {
        string name = "x[" + to_string(i) + "," + to_string(c) + "]";
        x[{i, c}] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        cnt++;
      }
    }
  }
  cout << "created " << cnt << " x variables" << endl;
}

inline void var_I(GRBModel *model, GRBVar *&I, vector<int> &conns) {
  I = new GRBVar[N];
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    string name = "I[" + to_string(i) + "]";
    I[i] = model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
  }
}

inline void var_Iij(GRBModel *model, mt2 &Iij, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsingChannel(i, c)) {
        string name = "Iij[" + to_string(i) + "," + to_string(c) + "]";
        Iij[{i, c}] =
            model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
      }
    }
  }
}

inline void var_z(GRBModel *model, mt2 &z, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsingChannel(i, c)) {
        string name = "z[" + to_string(i) + "," + to_string(c) + "]";
        z[{i, c}] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
      }
    }
  }
}

// -------------- constraints ------------------

inline void couple(GRBModel *model, mt3 &y, mt2 &x, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int b = 0; b < 4; ++b) {
      bool one = false, two = false;
      GRBLinExpr expr = 0, expr1 = 0;

      for (int m = 0; m < 12; ++m)
        if (canTransmitUsingBandwidth(i, b, m)) {
          one = true;
          assert(y.find({i, b, m}) != y.end());
          expr += y[{i, b, m}];
        }

      for (int c = 0; c < C_b[b].size(); ++c)
        if (canTransmitUsingChannel(i, C_b[b][c])) {
          two = true;
          assert(x.find({i, C_b[b][c]}) != x.end());
          expr1 += x[{i, C_b[b][c]}];
        }

      assert(one == two);
      if (!one) {
        continue;
      }

      string name = "couple[" + to_string(i) + "]";
      model->addConstr(expr <= expr1, name);
    }
  }
}

inline void unique(GRBModel *model, mt2 &x, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      expr += x[{i, c}];
    }

    string name = "unique[" + to_string(i) + "]";
    model->addConstr(expr <= 1, name);
  }
}

inline void ch_overlap(GRBModel *model, mt2 &z, mt2 &x, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c1 = 0; c1 < C; ++c1) {
      if (!canTransmitUsingChannel(i, c1)) continue;

      GRBLinExpr expr = 0;
      for (int c2 = 0; c2 < C; ++c2) {
        if (overlap[c1][c2] && canTransmitUsingChannel(i, c2))
          expr += x[{i, c2}];
      }

      string name = "over[" + to_string(i) + "," + to_string(c1) + "]";
      model->addConstr(expr == z[{i, c1}], name);
    }
  }
}

inline void interch(GRBModel *model, mt2 &z, mt2 &Iij, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      GRBLinExpr expr = 0;
      for (int _u = 0; _u < conns.size(); ++_u) {
        int u = conns[_u];
        if (u != i && canTransmitUsingChannel(u, c))
          expr += AFF[u][i] * z[{u, c}];
      }

      string name = "interch[" + to_string(i) + "," + to_string(c) + "]";
      model->addConstr(expr == Iij[{i, c}], name);
    }
  }
}

inline void bigG(GRBModel *model, GRBVar *I, mt2 &Iij, mt2 &x,
                 vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      string name = "bigG[" + to_string(i) + "," + to_string(c) + "]";
      model->addConstr(I[i] >= Iij[{i, c}] - BM[i] * (1 - x[{i, c}]), name);
    }
  }
}

inline void bigL(GRBModel *model, GRBVar *I, mt2 &Iij, mt2 &x,
                 vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      string name = "bigL[" + to_string(i) + "," + to_string(c) + "]";
      model->addConstr(I[i] <= Iij[{i, c}] + BM[i] * (1 - x[{i, c}]), name);
    }
  }
}

inline void sinr(GRBModel *model, GRBVar *I, mt2 &x, vector<int> &conns) {
  for (int _i = 0; _i < conns.size(); ++_i) {
    int i = conns[_i];
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      expr += (AFF[i][i] / B[i][cToBIdx(c)] - NOI) * x[{i, c}];
    }

    string name = "sinr[" + to_string(i) + "]";
    model->addConstr(I[i] <= expr, name);
  }
}

int main(int argc, char **argv) {
  freopen(argv[3], "r", stdin);
  printf("opening %s\n", argv[3]);
  read_data();
  set<int> scheduled;
  int off = N, OF = 0;

  double elapsed = 0.0;
  string solFileName = "sol" + string(argv[1]) + "_" + string(argv[2]) + ".sol";
  FILE *solFile = fopen(solFileName.c_str(), "a");
  while (off > 0) {
    auto start = std::chrono::steady_clock::now();
    OF++;

    vector<int> l_to_idx;
    for (int i = 0, idx = l_to_idx.size(); i < N; ++i)
      if (scheduled.find(i) == scheduled.end()) l_to_idx.push_back(i);

    GRBEnv env;
    GRBModel *model = new GRBModel(env);
    GRBVar *I, *t;
    mt2 x, Iij, z;
    mt3 y;
    // variables
    // printf("variables...\n");
    var_x(model, x, l_to_idx);
    var_z(model, z, l_to_idx);
    var_Iij(model, Iij, l_to_idx);
    var_I(model, I, l_to_idx);
    var_y(model, y, l_to_idx);

    // constraints
    model->update();
    // printf("constraints...\n");
    unique(model, x, l_to_idx);
    couple(model, y, x, l_to_idx);
    ch_overlap(model, z, x, l_to_idx);
    interch(model, z, Iij, l_to_idx);
    bigG(model, I, Iij, x, l_to_idx);
    bigL(model, I, Iij, x, l_to_idx);
    sinr(model, I, x, l_to_idx);
    // optimize

    model->set(GRB_DoubleParam_TimeLimit,
               min(600.0, max(0.0, 3600.0 - elapsed)));
    model->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
    // model->set(GRB_IntParam_LogToConsole, 0);
    model->set(GRB_DoubleParam_IntFeasTol, 1e-5);
    model->update();
    model->optimize();

    if (model->get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
      bool at_least_one = false;
      for (auto seila : x) {
        auto &[i, c] = seila.first;
        if (seila.second.get(GRB_DoubleAttr_X) == 1.0) {
          off--;
          scheduled.insert(i);
          at_least_one = true;
          fprintf(solFile, "%d %d %d\n", i, c, OF);
        }
      }

      // cout << "off is " << off << endl;
      assert(at_least_one);
    } else {
      if (model->get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
        model->computeIIS();
        model->write("why.ilp");
      }
      cout << "ops!" << endl;
      OF = -1;
      break;
    }

    auto finish = std::chrono::steady_clock::now();
    double elapsed_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(finish -
                                                                  start)
            .count();

    elapsed += elapsed_seconds;

    cout << "fooo " << elapsed_seconds << endl;
    if (elapsed > 3605) {
      cout << "time limit" << endl;
      break;
    }
  }

  fclose(solFile);

  FILE *obj = fopen("obj", "a");
  fprintf(obj, "%s\t%d\t%.3lf\n", argv[2], OF, elapsed);
  fclose(obj);
  return 0;
}
