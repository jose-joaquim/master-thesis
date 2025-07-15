#include "basic.h"
#include "gurobi_c++.h"
#include "gurobi_c.h"

using namespace std;
using mt3 = map<tuple<int, int, int>, GRBVar>;

bool canTransmitUsing(int i, int c) {
  return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]) ||
         approximatelyEqual(DR[11][cToBIdx(c)], GMM[i]);
}

void build_subproblem(GRBModel &model, mt3 &x, vector<GRBVar> &y,
                      vector<double> pi) {
  // variables
  printf("variables...\n");

  for (int i = 0; i < N; ++i)
    y.push_back(model.addVar(0.0, 1.0, pi[i], GRB_BINARY, "y" + to_string(i)));

  for (int i = 0; i < N; ++i)
    for (int c = 0; c < C; ++c)
      for (int m = 0; m < M; ++m)
        if (canTransmitUsing(i, c))
          x[{i, c, m}] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                      "x" + to_string(i) + "," + to_string(c) +
                                          "," + to_string(m));

  // constraints
  model.update();

  printf("constraints...\n");
  for (int i = 0; i < N; ++i) {
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c)
      if (canTransmitUsing(i, c))
        for (int m = 0; m < M; ++m)
          expr += x[{i, c, m}];

    model.addConstr(expr == y[i], "y==sum(x)_" + to_string(i));
  }

  for (int i = 0; i < N; ++i)
    for (int c = 0; c < C; ++c)
      if (canTransmitUsing(i, c))
        for (int m = 0; m < M; ++m) {
          GRBLinExpr noise = -BM[i] * (1 - x[{i, c, m}]);

          for (int u = 0; u < N; ++u)
            for (int c_ = 0; c_ < C; ++c_)
              if (canTransmitUsing(u, c_))
                if (u != i && overlap[c][c_]) {
                  for (int m_ = 0; m_ < M; ++m_)
                    noise += AFF[u][i] * x[{u, c_, m_}];
                }

          double sinr_at_channel = AFF[i][i] / B[i][cToBIdx(c)] - NOI;
          model.addConstr(x[{i, c, m}] * sinr_at_channel >= noise,
                          "sinr>=noise_" + to_string(i) + "_" + to_string(c) +
                              "_" + to_string(m));
        }

  model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
}

int main(int argc, char *argv[]) {
  freopen(argv[1], "r", stdin);
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
      masterConstrs.emplace_back(
          master.addConstr(vars[l] >= 1.0, "link" + to_string(l)));

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
           << relax.get(GRB_DoubleAttr_ObjVal) << endl;
      vector<double> pi;
      for (int l = 0; l < N; ++l)
        pi.emplace_back(relax.getConstr(l).get(GRB_DoubleAttr_Pi));

      GRBModel sub(env);
      vector<GRBVar> y;
      mt3 x;
      build_subproblem(sub, x, y, pi);

      sub.write("sub" + to_string(loops) + ".lp");
      sub.set(GRB_IntParam_LogToConsole, 0);
      sub.optimize();

      if (sub.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
        cout << "sub solved! objective: " << sub.get(GRB_DoubleAttr_ObjVal)
             << endl;

      map<int, int> solution;
      for (int l = 0; l < N; ++l)
        solution[l] = static_cast<int>(y[l].get(GRB_DoubleAttr_X) + .5);

      map<int, GRBColumn> col;
      for (int l = 0; l < N; ++l)
        if (solution[l] > 0)
          col[l].addTerm(solution[l], masterConstrs[l]);

      for (auto &e : col)
        vars.emplace_back(
            master.addVar(0.0, 1.0, 1.0, GRB_BINARY, e.second,
                          "x" + to_string(loops) + "_" + to_string(e.first)));

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

    FILE *obj = fopen("obj", "a");
    fprintf(obj, "%lf\n", master.get(GRB_DoubleAttr_ObjVal));
    fclose(obj);
    cout << master.get(GRB_DoubleAttr_ObjVal) << endl;
  } catch (GRBException e) {
    cout << "Error code = " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  } catch (...) {
    cout << "Exception during optimization" << endl;
  }

  return 0;
}
