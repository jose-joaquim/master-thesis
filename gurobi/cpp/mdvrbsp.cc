#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <tuple>

#include "gurobi_c++.h"
#include "header.h"

using namespace std;
using mt3 = map<tuple<int, int, int>, GRBVar>;

bool canTransmitUsing(int i, int c) {
  return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]) ||
         approximatelyEqual(DR[11][cToBIdx(c)], GMM[i]);
}

// -------------- variables --------------------
inline void var_x(GRBModel *model, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsing(i, c)) {
        for (int t = 0; t < T; ++t) {
          string name = "x[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
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
          string name = "Iij[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
          Iij[{i, c, t}] =
              model->addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, name);
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
          string name = "z[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
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
  cout << "unique" << endl;
  for (int i = 0; i < N; ++i) {
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c)) continue;

      for (int t = 0; t < T; ++t) expr += x[{i, c, t}];
    }

    string name = "unique[" + to_string(i) + "]";
    model->addConstr(expr == 1, name);
  }
}

inline void waste(GRBModel *model, mt3 &x, GRBVar *t) {
  cout << "waste" << endl;
  for (int i = 0; i < N; ++i) {
    for (int _t = 0; _t < T; ++_t) {
      GRBLinExpr expr = 0;
      for (int c = 0; c < C; ++c)
        if (canTransmitUsing(i, c)) expr += x[{i, c, _t}];

      string name = "waste[" + to_string(i) + "," + to_string(_t) + "]";
      model->addConstr(expr <= t[_t], name);
    }
  }
}

inline void ch_overlap(GRBModel *model, mt3 &z, mt3 &x) {
  cout << "overlap" << endl;
  for (int i = 0; i < N; ++i) {
    for (int c1 = 0; c1 < C; ++c1) {
      if (!canTransmitUsing(i, c1)) continue;

      for (int t = 0; t < T; ++t) {
        GRBLinExpr expr = 0;
        for (int c2 = 0; c2 < C; ++c2) {
          if (overlap[c1][c2] && canTransmitUsing(i, c2)) expr += x[{i, c2, t}];
        }

        string name = "over[" + to_string(i) + "," + to_string(c1) + "," +
                      to_string(t) + "]";
        model->addConstr(expr == z[{i, c1, t}], name);
      }
    }
  }
}

inline void interch(GRBModel *model, mt3 &z, mt3 &Iij) {
  cout << "interch" << endl;
  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < T; ++t) {
      for (int c = 0; c < C; ++c) {
        if (!canTransmitUsing(i, c)) continue;

        GRBLinExpr expr = 0;
        for (int u = 0; u < N; ++u) {
          if (u != i && canTransmitUsing(u, c))
            expr += AFF[u][i] * z[{u, c, t}];
        }

        string name = "interch[" + to_string(i) + "," + to_string(c) + "," +
                      to_string(t) + "]";
        model->addConstr(expr == Iij[{i, c, t}], name);
      }
    }
  }
}

inline void bigG(GRBModel *model, GRBVar *I, mt3 &Iij, mt3 &x) {
  cout << "bigG" << endl;
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c)) continue;

      for (int t = 0; t < T; ++t) {
        string name = "bigG[" + to_string(i) + "," + to_string(c) + "," +
                      to_string(t) + "]";
        model->addConstr(I[i] >= Iij[{i, c, t}] - BM[i] * (1 - x[{i, c, t}]),
                         name);
      }
    }
  }
}

inline void bigL(GRBModel *model, GRBVar *I, mt3 &Iij, mt3 &x) {
  cout << "bigL" << endl;
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c)) continue;

      for (int t = 0; t < T; ++t) {
        string name = "bigL[" + to_string(i) + "," + to_string(c) + "," +
                      to_string(t) + "]";
        model->addConstr(I[i] <= Iij[{i, c, t}] + BM[i] * (1 - x[{i, c, t}]),
                         name);
      }
    }
  }
}

inline void sinr(GRBModel *model, GRBVar *I, mt3 &x) {
  cout << "sinr" << endl;
  for (int i = 0; i < N; ++i) {
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsing(i, c)) continue;

      for (int t = 0; t < T; ++t) {
        expr += (AFF[i][i] / B[i][cToBIdx(c)] - NOI) * x[{i, c, t}];
      }
    }

    string name = "sinr[" + to_string(i) + "]";
    model->addConstr(I[i] <= expr, name);
  }
}

inline void pairwise(GRBModel *model, mt3 &x, mt3 &z) {
  cout << "pairwise" << endl;
  using ptv = pair<tuple<int, int, int>, GRBVar>;
  vector<vector<ptv>> var_l(N, vector<ptv>());

  for (const auto &[key_x, var_x] : x) {
    const auto &[l_x, c_x, t_x] = key_x;
    var_l[l_x].push_back(make_pair(key_x, var_x));
  }

  for (int i = 0; i < N; ++i) {
    for (const auto &[key_i, var_i] : var_l[i]) {
      const auto &[_, c_i, t_i] = key_i;
      for (int j = 0; j < N; ++j) {
        if (i == j) continue;

        int bw_idx = cToBIdx(c_i);
        double aff = B[i][bw_idx] * AFF[j][i] / AFF[i][i];
        if (definitelyGreaterThan(aff, 1.0)) {
          GRBLinExpr expr = var_i;
          bool any = false;

          for (const auto &[key_j, var_j] : var_l[j]) {
            const auto &[l_j, c_j, t_j] = key_j;

            int bw_idx_2 = cToBIdx(c_j);
            if (!overlap[c_i][c_j] || t_i != t_j || bw_idx != bw_idx_2)
              continue;

            any = true;
            expr += var_j;
          }

          if (any) model->addConstr(expr <= 1);
        }
      }
    }
  }

  // for (const auto &[key_x, var_x] : x) {
  //   const auto &[l_x, c_x, t_x] = key_x;
  //
  //   for (int k = 0; k < N; ++k) {
  //     if (l_x == k) continue;
  //
  //     int bw_idx = cToBIdx(c_x);
  //     double aff = B[l_x][bw_idx] * AFF[k][l_x] / AFF[l_x][l_x];
  //     if (definitelyGreaterThan(aff, 1.0)) {
  //       GRBLinExpr expr = var_x;
  //       bool any = false;
  //
  //       for (const auto &[key_y, var_y] : x) {
  //         const auto &[l_y, c_y, t_y] = key_y;
  //
  //         int bw_idx_2 = cToBIdx(c_y);
  //         if (t_y > l_x || l_x == l_y || l_y != k || !overlap[c_x][c_y] ||
  //             t_x != t_y || bw_idx != bw_idx_2)
  //           continue;
  //
  //         any = true;
  //         expr += var_y;
  //       }
  //
  //       if (any) model->addConstr(expr <= 1);
  //     }
  //   }
  // }
}

void fix_variables(GRBModel *model, mt3 &x, string file_name) {
  char name[100];
  double value = 0;
  FILE *file = fopen(file_name.c_str(), "r");
  fgets(name, 100, file);
  while (fscanf(file, "%s %lf\n", name, &value) != EOF) {
    if (name[0] != 'x') continue;
    GRBVar x = model->getVarByName(string(name));
    model->addConstr(x >= value);
    // printf("li %s %lf\n", name, value);
  }

  fclose(file);
}

class LogCallback : public GRBCallback {
 public:
  string path;
  vector<string> lines;
  double last_call = 0;

  LogCallback(string path) : path(path) {}
  ~LogCallback() {
    ofstream file{path};

    if (!file.is_open()) {
      cout << "nao consegui abrir o arquivo callback" << endl;
      exit(10);
    }

    for (string line : lines) file << line;
    file.close();
  }

 protected:
  void callback() {
    try {
      if (where == GRB_CB_MIP) {
        double runtime = getDoubleInfo(GRB_CB_RUNTIME);

        if (approximatelyEqual(last_call, 0.0) ||
            definitelyGreaterThan(runtime - last_call, 5)) {
          last_call = runtime;
          double objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
          if (objbst == GRB_INFINITY) objbst = 0.0;
          double objbnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
          if (objbnd == GRB_INFINITY) objbst = 0.0;

          double nodecnt = getDoubleInfo(GRB_CB_MIP_NODCNT);
          int solcnt = getIntInfo(GRB_CB_MIP_SOLCNT);

          char line[100];
          snprintf(line, 100, "%.3lf %.3lf %.3lf %.3lf %d\n", runtime, nodecnt,
                   objbst, objbnd, solcnt);
          lines.push_back(string(line));
        }
      } /* else if (where == GRB_CB_MIPNODE) {
        int status = getIntInfo(GRB_CB_MIPNODE_STATUS);

        if (status == GRB_OPTIMAL) {
          double runtime = getDoubleInfo(GRB_CB_RUNTIME);
          double objbst = getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
          int solcnt = getIntInfo(GRB_CB_MIPNODE_SOLCNT);

          char line[100];
          snprintf(line, 100, "%.3lf %.3lf %.3lf %d\n", runtime, objbst, objbst,
                   solcnt);
          lines.push_back(string(line));
          }
          }*/
    } catch (GRBException e) {
      cout << "Error number: " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    } catch (...) {
      cout << "Error during callback" << endl;
    }
  }
};

// ---------------------------------------------

double run(char **argv, bool pair) {
  auto start = high_resolution_clock::now();
  GRBEnv env;
  GRBModel *model = new GRBModel(env);
  GRBVar *I, *t;
  mt3 x, Iij, z;
  // variables
  printf("variables...\n");
  var_x(model, x);
  var_z(model, z);
  var_Iij(model, Iij);
  var_I(model, I);
  var_t(model, t);

  auto stop = high_resolution_clock::now();
  duration<double> ms_double = stop - start;
  cout << ms_double.count() << endl;

  // constraints
  model->update();
  printf("constraints...\n");
  symmetry1(model, t);
  symmetry2(model, x);
  unique(model, x);
  waste(model, x, t);
  ch_overlap(model, z, x);
  interch(model, z, Iij);
  bigG(model, I, Iij, x);
  bigL(model, I, Iij, x);
  sinr(model, I, x);
  if (pair) pairwise(model, x, z);

  stop = high_resolution_clock::now();
  ms_double = stop - start;
  cout << ms_double.count() << endl;
  // optimize
  model->update();

  model->set(GRB_IntParam_LogToConsole, 0);
  model->set(GRB_DoubleParam_IntFeasTol, 1e-7);

  string log_cb =
      "./callback_log" + string(argv[1]) + "_" + string(argv[2]) + ".txt";
  LogCallback cb = LogCallback(log_cb);

  model->setCallback(&cb);
  model->optimize();

  if (model->get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
    string res = "result" + string(argv[1]);
    FILE *result = fopen(res.c_str(), "a");
    double objVal = -1;
    if (model->get(GRB_IntAttr_Status) != GRB_TIME_LIMIT)
      objVal = model->get(GRB_DoubleAttr_ObjVal);

    double dualObj = model->get(GRB_DoubleAttr_ObjBoundC);

    if (approximatelyEqual(objVal, GRB_INFINITY * 1.0)) objVal = -1.0;

    if (approximatelyEqual(dualObj, GRB_INFINITY * 1.0)) dualObj = -1.0;

    model->write(argv[4]);

    fprintf(result, "%s\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\n", argv[2],
            objVal, dualObj, model->get(GRB_DoubleAttr_MIPGap),
            model->get(GRB_IntAttr_NumVars), model->get(GRB_IntAttr_NumConstrs),
            model->get(GRB_DoubleAttr_DNumNZs),
            model->get(GRB_DoubleAttr_NodeCount),
            model->get(GRB_DoubleAttr_Runtime));
    fclose(result);

    string res_sol = "sol" + string(argv[1]) + "_" + string(argv[2]) + ".sol";
    FILE *sol = fopen(res_sol.c_str(), "w");

    for (auto &var : x) {
      const auto &[i, c, t] = var.first;
      if (approximatelyEqual(var.second.get(GRB_DoubleAttr_X), 1.0))
        fprintf(sol, "%d %d %d %lf\n", i, c, t,
                var.second.get(GRB_DoubleAttr_X));
    }

    fclose(sol);
    return objVal;
  } else {
    model->set(GRB_IntParam_IISMethod, 1);
    model->computeIIS();
    model->write("xi.ilp");
    cout << "not optimal!" << endl;
    return -1;
  }
}

int main(int argc, char **argv) {
  freopen(argv[3], "r", stdin);
  read_data();

  double obj = run(argv, true);
  return 0;
}
