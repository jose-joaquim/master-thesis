#! /usr/bin/env python3
# With gurobi syntatic sugars, with Big-M
# This model has a dictionary to select the frequency spectrum
# that is allowed to be used

import gurobipy as gp
from gurobipy import GRB

# connections, time-slots, big-M, affectance, channels dictionary
def defineModel(N, NTS, BM, AFF, AUXNC, B, NOI, OVER, cToBIdx, DR, GMM):
    m = gp.Model("mdvrbsp reduced channels")

    def variables():
        t = m.addVars(NTS, vtype=GRB.BINARY, name="t")
        x = m.addVars(N, len(AUXNC), NTS, vtype=GRB.BINARY, name="x")
        I_ = m.addVars(N, lb=0.0, vtype=GRB.CONTINUOUS, name="I")
        Iij = m.addVars(N, len(AUXNC), NTS, lb=0.0, vtype=GRB.CONTINUOUS, name="Iij")
        z = m.addVars(N, len(AUXNC), NTS, vtype=GRB.BINARY, name="z")

        return t, x, I_, Iij, z

    def constraints(t, x, I_, Iij, z):
        # 0
        m.addConstrs(
            (
                x[i, c, t] == 0
                for c in AUXNC
                for i in range(N)
                for t in range(NTS)
                if DR[11][cToBIdx(c)] < GMM[i]
            ),
            "cant",
        )

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
                gp.quicksum(AFF[u][i] * z[u, c, t_] for u in range(N) if u != i)
                == Iij[i, c, t_]
                for t_ in range(NTS)
                for c in AUXNC
                for i in range(N)
            ),
            "interch",
        )

        # 6.1
        m.addConstrs(
            (
                I_[i] >= Iij[i, c, t_] - BM[i] * (1 - x[i, c, t_])
                for t_ in range(NTS)
                for c in AUXNC
                for i in range(N)
            ),
            "bigG_",
        )

        # 6.2
        m.addConstrs(
            (
                I_[i] <= Iij[i, c, t_] + BM[i] * (1 - x[i, c, t_])
                for t_ in range(NTS)
                for c in AUXNC
                for i in range(N)
            ),
            "bigL_",
        )

        # 7
        m.addConstrs(
            (
                I_[i]
                <= gp.quicksum(
                    (AFF[i][i] / B[i][cToBIdx(c[1])] - NOI) * x[i, c[0], t_]
                    for t_ in range(NTS)
                    for c in AUXNC.items()
                )
                for i in range(N)
            ),
            "sinr",
        )

    def objective(t):
        m.setObjective(t.sum(), GRB.MINIMIZE)

    t, x, I_, Iij, z = variables()
    constraints(t, x, I_, Iij, z)
    objective(t)
    return m, x, I_
