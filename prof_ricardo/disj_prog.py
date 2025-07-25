#! /usr/bin/env python3
# With gurobi syntatic sugars, with Big-M
# This model has a dictionary to select the frequency spectrum
# that is allowed to be used

import gurobipy as gp
import math
from gurobipy import GRB


def defineModelNonLinear(
    N,  # nConnections
    B,  # set of channel bandwidths (MHz) [20,40,80,160]
    NC,  # set of subdivisions of the electromagnetic spectrum into communication channels
    MCS,  # Modulation and Coding Scheme
    OVER,  # overlap
    DR,  # datarates
    NS,  # noise
    AFF,  # affectance
    SINR,  # signal to interference plus noise ratio
    BM,
):  # big M
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

    m = gp.Model("nonlinear")

    def variables():
        x = m.addVars(N, NC, vtype=GRB.BINARY, name="x")
        y = m.addVars(N, B, MCS, vtype=GRB.BINARY, name="y")
        # y = m.addVars(N, B, MCS, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")

        return x, y

    def constraints(x, y):
        # Constraint 1
        m.addConstrs(
            (gp.quicksum(x[i, c] for c in range(NC)) <= 1 for i in range(N)), "unique"
        )
        # Constraint 2
        m.addConstrs(
            gp.quicksum(
                y[i, b, m] * ((AFF[i][i] / SINR[m][b]) - NS)
                for b in range(B)
                for m in range(MCS)
            )
            >= gp.quicksum(
                AFF[u][i] * x[u, _c] * x[i, c]
                for u in range(N)
                for c in range(NC)
                for _c in range(NC)
                if OVER[c][_c] == 1 and i != u
            )
            for i in range(N)
        )
        # Constraint 4
        m.addConstrs(
            gp.quicksum(y[i, b, m] for m in range(MCS))
            == gp.quicksum(x[i, c] for c in C_b[b])
            for b in range(B)
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

    x, y = variables()
    constraints(x, y)
    objective(y)
    return m, x, y


def defineModelDisj(
    N,  # nConnections
    B,  # set of channel bandwidths (MHz) [20,40,80,160]
    NC,  # set of subdivisions of the electromagnetic spectrum into communication channels
    MCS,  # Modulation and Coding Scheme
    OVER,  # overlap
    DR,  # datarates
    NS,  # noise
    AFF,  # affectance
    SINR,  # signal to interference plus noise ratio
    BM,
):  # big M
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

    m = gp.Model("disj")

    def variables():
        x = m.addVars(N, NC, vtype=GRB.BINARY, name="x")
        y = m.addVars(N, B, MCS, vtype=GRB.BINARY, name="y")
        # y = m.addVars(N, B, MCS, lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS, name="y")

        return x, y

    def constraints(x, y):
        # Constraint 1
        m.addConstrs(
            (gp.quicksum(x[i, c] for c in range(NC)) <= 1 for i in range(N)), "unique"
        )
        # Constraint 2
        # BigM = [max(((AFF[i][i] / SINR[m][b]) - NS) for b in range(B) for m in range(MCS)) for i in range(N)]
        BigM = [max(AFF[u][i] for u in range(N)) for i in range(N)]

        m.addConstrs(
            gp.quicksum(
                y[i, b, m] * ((AFF[i][i] / SINR[m][b]) - NS)
                for b in range(B)
                for m in range(MCS)
            )
            >= gp.quicksum(
                AFF[u][i] * x[u, _c]
                for u in range(N)
                for _c in range(NC)
                if OVER[c][_c] == 1 and i != u
            )
            - BigM[i] * (1 - x[i, c])
            for i in range(N)
            for c in range(NC)
        )
        # Constraint 4
        m.addConstrs(
            gp.quicksum(y[i, b, m] for m in range(MCS))
            == gp.quicksum(x[i, c] for c in C_b[b])
            for b in range(B)
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

    x, y = variables()
    constraints(x, y)
    objective(y)
    return m, x, y
