#! /usr/bin/env python3

import gurobipy as gp
from gurobipy import GRB


def defineModel(N, NTS, AFF, AUXNC, B, NOI, OVER, cToBIdx):
    m = gp.Model("mdvrbsp quadratic variables formulation")

    def variables():
        t = m.addVars(NTS, vtype=GRB.BINARY, name="t")
        x = m.addVars(N, len(AUXNC), NTS, vtype=GRB.BINARY, name="x")
        I_ = m.addVars(N, lb=0.0, name="I_")
        z = m.addVars(N, len(AUXNC), NTS, vtype=GRB.BINARY, name="z")
        w = m.addVars(N, N, name="w")

        return t, x, I_, z, w

    def constraints(t, x, I_, z, w):
        # 0
        m.addConstr(t.sum() >= 1, "force")

        # 1
        m.addConstrs((t[i + 1] <= t[i] for i in range(NTS - 1)), "ts")

        # 2
        m.addConstrs(
            (
                gp.quicksum(x[i, c, t_] for c in AUXNC for t_ in range(NTS)) == 1
                for i in range(N)
            ),
            "unique",
        )

        # 3
        m.addConstrs(
            (
                gp.quicksum(x[i, c, t_] for c in AUXNC) <= t[t_]
                for t_ in range(NTS)
                for i in range(N)
            ),
            "waste",
        )

        # 4
        m.addConstrs(
            (
                gp.quicksum(
                    x[i, c2[0], t_] for c2 in AUXNC.items() if OVER[c1[1]][c2[1]] == 1
                )
                == z[i, c1[0], t_]
                for t_ in range(NTS)
                for c1 in AUXNC.items()
                for i in range(N)
            ),
            "overlap",
        )

        # 5
        m.addConstrs(
            (
                w[i, u] >= x[i, c, t_] + z[u, c, t_] - 1
                for c in AUXNC
                for t_ in range(NTS)
                for u in range(N)
                for i in range(N)
                if u != i
            ),
            "same",
        )

        # 6
        m.addConstrs(
            (
                gp.quicksum(AFF[u][i] * w[u, i] for u in range(N) if u != i) == I_[i]
                for i in range(N)
            ),
            "inter",
        )

        # 7
        m.addConstrs(
            (
                I_[i]
                <= gp.quicksum(
                    (AFF[i][i] / B[i][cToBIdx(c[1])] * x[i, c[0], t_]) - NOI
                    for t_ in range(NTS)
                    for c in AUXNC.items()
                )
                for i in range(N)
            ),
            "sinr",
        )

    def objective(t):
        m.setObjective(t.sum(), GRB.MINIMIZE)

    t, x, I_, z, w = variables()
    constraints(t, x, I_, z, w)
    objective(t)
    return m, x, I_
