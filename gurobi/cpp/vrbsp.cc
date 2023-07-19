#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <tuple>

#include "gurobi_c++.h"
#include "header.h"

using namespace std;
using mt3 = map<tuple<int, int, int>, GRBVar>;

bool canTransmitUsingChannel(int i, int c) {
  return definitelyGreaterThan(DR[11][cToBIdx(c)], GMM[i]);
}

bool canTransmitUsingBandwidth(int i, int b, int m) {
  return definitelyGreaterThan(DR[m][b], GMM[i]);
}

// -------------- variables --------------------
inline void var_y(GRBModel *model, mt3 y) {
  for (int i = 0; i < N; ++i) {
    for (int b = 0; b < 4; ++b) {
      for (int m = 0; m < 12; ++m) {
        string name =
            "y[" + to_string(i) + "," + to_string(b) + "," + to_string(m) + "]";
        y[{i, b, m}] = model->addVar(0.0, 1.0, DR[m][b], GRB_BINARY, name);
      }
    }
  }
}

inline void var_x(GRBModel *model, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (canTransmitUsingChannel(i, c)) {
        for (int t = 0; t < T; ++t) {
          string name = "x[" + to_string(i) + "," + to_string(c) + "," +
                        to_string(t) + "]";
          x[{i, c, t}] = model->addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
        }
      }
    }
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
      if (canTransmitUsingChannel(i, c)) {
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
      if (canTransmitUsingChannel(i, c)) {
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

inline void unique(GRBModel *model, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      for (int t = 0; t < T; ++t) expr += x[{i, c, t}];
    }

    string name = "unique[" + to_string(i) + "]";
    model->addConstr(expr == 1, name);
  }
}

inline void ch_overlap(GRBModel *model, mt3 &z, mt3 &x) {
  for (int i = 0; i < N; ++i) {
    for (int c1 = 0; c1 < C; ++c1) {
      if (!canTransmitUsingChannel(i, c1)) continue;

      for (int t = 0; t < T; ++t) {
        GRBLinExpr expr = 0;
        for (int c2 = 0; c2 < C; ++c2) {
          if (overlap[c1][c2] && canTransmitUsingChannel(i, c2))
            expr += x[{i, c2, t}];
        }

        string name = "over[" + to_string(i) + "," + to_string(c1) + "," +
                      to_string(t) + "]";
        model->addConstr(expr == z[{i, c1, t}], name);
      }
    }
  }
}

inline void interch(GRBModel *model, mt3 &z, mt3 &Iij) {
  for (int i = 0; i < N; ++i) {
    for (int t = 0; t < T; ++t) {
      for (int c = 0; c < C; ++c) {
        if (!canTransmitUsingChannel(i, c)) continue;

        GRBLinExpr expr = 0;
        for (int u = 0; u < N; ++u) {
          if (u != i && canTransmitUsingChannel(u, c))
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
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

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
  for (int i = 0; i < N; ++i) {
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

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
  for (int i = 0; i < N; ++i) {
    GRBLinExpr expr = 0;
    for (int c = 0; c < C; ++c) {
      if (!canTransmitUsingChannel(i, c)) continue;

      for (int t = 0; t < T; ++t) {
        expr += (AFF[i][i] / B[i][cToBIdx(c)] - NOI) * x[{i, c, t}];
      }
    }

    string name = "sinr[" + to_string(i) + "]";
    model->addConstr(I[i] <= expr, name);
  }
}

int main() { return 0; }
