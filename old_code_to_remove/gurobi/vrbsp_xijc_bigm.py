#! /usr/bin/env python3
# With gurobi syntatic sugars, with Big-M
# This model has a dictionary to select the frequency spectrum
# that is allowed to be used

import gurobipy as gp
import math
from gurobipy import GRB


def defineModel(N, B, NC, MCS, OVER, DR, NS, AFF, SINR, BM, cToB):
    m = gp.Model("vrbsp xijcm bigm")

    def variables():
        x = m.addVars(N, NC, MCS, vtype=GRB.BINARY, name="x")
        z = m.addVars(N, NC, vtype=GRB.BINARY, name="z")
        I_ = m.addVars(N, lb=0.0, vtype=GRB.CONTINUOUS, name="I")
        Iij = m.addVars(N, NC, lb=0.0, vtype=GRB.CONTINUOUS, name="Iij")

        return x, z, I_, Iij

    def constraints(x, z, I_, Iij):
        # Constraint 1
        m.addConstrs(
            (
                gp.quicksum(x[i, c, m] for c in range(NC) for m in range(MCS)) <= 1.0
                for i in range(N)
            ),
            "unique",
        )

        # Constraint 2
        m.addConstrs(
            (
                gp.quicksum(
                    x[i, c2, m]
                    for c2 in range(NC)
                    if OVER[c1][c2] == 1
                    for m in range(MCS)
                )
                == z[i, c1]
                for c1 in range(NC)
                for i in range(N)
            ),
            "overlap",
        )

        # Constraint 3
        m.addConstrs(
            (
                gp.quicksum(AFF[u][i] * z[u, c] for u in range(N) if i != u)
                == Iij[i, c]
                for c in range(NC)
                for i in range(N)
            ),
            "aff",
        )

        # Constraint 4
        m.addConstrs(
            (
                I_[i] <= Iij[i, c] + BM[i] * (1 - x[i, c, m])
                for c in range(NC)
                for m in range(MCS)
                for i in range(N)
            ),
            "bigL",
        )

        # Constraint 5
        m.addConstrs(
            (
                I_[i] >= Iij[i, c] - BM[i] * (1 - x[i, c, m])
                for c in range(NC)
                for m in range(MCS)
                for i in range(N)
            ),
            "bigG",
        )

        # Constraint 6
        m.addConstrs(
            (
                I_[i]
                <= gp.quicksum(
                    x[i, c, m] * ((AFF[i][i] / SINR[m][cToB(c)]) - NS)
                    for c in range(NC)
                    for m in range(MCS)
                )
                for i in range(N)
            ),
            "sinr",
        )

    def objective(x):
        m.setObjective(
            sense=GRB.MAXIMIZE,
            expr=(
                gp.quicksum(
                    DR[m][cToB(c)] * x[i, c, m]
                    for c in range(NC)
                    for m in range(MCS)
                    for i in range(N)
                )
            ),
        )

    x, z, I_, Iij = variables()
    constraints(x, z, I_, Iij)
    objective(x)
    return m, x, I_
