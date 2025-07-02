#! /usr/bin/env python3
# With gurobi syntatic sugars, with Big-M
# This model has a dictionary to select the frequency spectrum
# that is allowed to be used

import gurobipy as gp
import math
from gurobipy import GRB


def defineModel(N, B, NC, MCS, OVER, DR, NS, AFF, SINR, BM):
    C_b = [
        [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            11,
            12,
            13,
            14,
            15,
            16,
            17,
            18,
            19,
            20,
            21,
            22,
            23,
            24,
        ],
        [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
        [37, 38, 39, 40, 41, 42],
        [43, 44],
    ]

    m = gp.Model("vrbsp ybm bigm")

    def variables():
        x = m.addVars(N, NC, vtype=GRB.BINARY, name="x")
        y = m.addVars(N, B, MCS, vtype=GRB.BINARY, name="y")
        #y = m.addVars(N, B, MCS, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")
        z = m.addVars(N, NC, vtype=GRB.BINARY, name="z")
        I_ = m.addVars(N, lb=0.0, vtype=GRB.CONTINUOUS, name="I")
        Iij = m.addVars(N, NC, lb=0.0, vtype=GRB.CONTINUOUS, name="Iij")

        return x, y, z, I_, Iij

    def constraints(x, y, z, I_, Iij):
        # Constraint 1
        m.addConstrs(
            (gp.quicksum(x[i, c] for c in range(NC)) <= 1.0 for i in range(N)), "unique"
        )

        # Constraint 2
        m.addConstrs(
            gp.quicksum(y[i, b, m] for m in range(MCS))
            <= gp.quicksum(x[i, c] for c in C_b[b])
            for b in range(B)
            for i in range(N)
        )

        # Constraint 3
        m.addConstrs(
            gp.quicksum(x[i, c2] for c2 in range(NC) if OVER[c1][c2] == 1) == z[i, c1]
            for c1 in range(NC)
            for i in range(N)
        )

        # Constraint 4
        m.addConstrs(
            gp.quicksum(AFF[u][i] * z[u, c] for u in range(N) if i != u) == Iij[i, c]
            for c in range(NC)
            for i in range(N)
        )

        # Constraint 5
        m.addConstrs(
            I_[i] <= Iij[i, c] + BM[i] * (1 - x[i, c])
            for c in range(NC)
            for i in range(N)
        )

        # Constraint 6
        m.addConstrs(
            I_[i] >= Iij[i, c] - BM[i] * (1 - x[i, c])
            for c in range(NC)
            for i in range(N)
        )

        # Constraint 7
        m.addConstrs(
            I_[i]
            <= gp.quicksum(
                y[i, b, m] * ((AFF[i][i] / SINR[m][b]) - NS)
                for b in range(B)
                for m in range(MCS)
            )
            for i in range(N)
        )

    def objective(y):
        m.setObjective(
            sense=GRB.MAXIMIZE,
            expr=(
                gp.quicksum(
                    DR[m][b] * y[i, b, m]
                    for b in range(B)
                    for m in range(MCS)
                    for i in range(N)
                )
            ),
        )

    x, y, z, I_, Iij = variables()
    constraints(x, y, z, I_, Iij)
    objective(y)
    return m, x, y, I_
