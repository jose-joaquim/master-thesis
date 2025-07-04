import gurobipy as gp
import math
from gurobipy import GRB


def defineAlternative1(
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
):
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

    def c_to_b(c):
        if c < 25:
            return 0
        elif c < 37:
            return 1
        elif c < 43:
            return 2
        return 3

    BigM = [max(AFF[u][i] for u in range(N)) for i in range(N)]
    m = gp.Model("alternative1")

    def variables():
        x = m.addVars(N, NC, vtype=GRB.BINARY, name="x")
        y = m.addVars(N, B, MCS, vtype=GRB.BINARY, name="y")
        w = m.addVars(N, N, NC, NC, vtype=GRB.BINARY, name="w")
        z = m.addVars(N, B, NC, MCS, vtype=GRB.BINARY, name="z")

        return x, y, w, z

    def constraints(x, y):
        m.addConstrs(
            (gp.quicksum(x[i, c] for c in range(NC)) <= 1.0 for i in range(N)), "unique"
        )

        m.addConstrs(
            (
                gp.quicksum(y[i, b, m] for m in range(MCS))
                <= gp.quicksum(x[i, c] for c in C_b[b])
                for b in range(B)
                for i in range(N)
            ),
            "couple",
        )

        m.addConstrs(
            (
                gp.quicksum(
                    (AFF[i][i] / SINR[m][c_to_b(c)]) * y[i, c_to_b(c), m]
                    for m in range(MCS)
                )
                >= gp.quicksum(
                    AFF[u][i] * x[u, c_]
                    for c_ in range(NC)
                    if OVER[c][c_]
                    for u in range(N)
                    if u != i
                )
                - BigM[i] * (1 - x[i, c])
                for i in range(N)
                for c in range(NC)
            ),
            "sinr",
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

    x, y, w, z = variables()
    constraints(x, y)
    objective(y)
    return m, x, y
