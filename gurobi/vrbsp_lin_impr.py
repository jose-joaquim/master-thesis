#! /usr/bin/env python3

import sys

import gurobipy as gp
import math
from gurobipy import GRB

overlap = []
nChan = 45
count_inst = 0
rst_headers = [
    "ObjVal",
    "ObjBoundC",
    "MIPGap",
    "NumVars",
    "NumConstrs",
    "IterCount",
    "BarIterCount",
    "NodeCount",
    "Runtime",
]


def convertTableToMW(SINR, SINR_):
    for i in range(len(SINR_)):
        for j in range(len(SINR_[i])):
            if SINR[i][j] != 0.0:
                b = SINR[i][j] / 10.0
                result = math.pow(10.0, b)

                SINR_[i][j] = result
            else:
                SINR_[i][j] = 0.0


def distance(a, b, c, d):
    return math.hypot((a - c), (b - d))


def distanceAndInterference(
    senders, receivers, interMatrix, distMatrix, powerSender, nConn, alfa,
):
    for i in range(nConn):
        X_si = receivers[i][0]
        Y_si = receivers[i][1]

        for j in range(nConn):
            X_rj = senders[j][0]
            Y_rj = senders[j][1]

            dist = distance(X_si, Y_si, X_rj, Y_rj)
            distMatrix[i][j] = dist

            if i == j:
                interMatrix[i][j] = 0.0
            else:
                value = powerSender / math.pow(dist, alfa) if dist != 0.0 else 1e9
                interMatrix[i][j] = value


def cToB(c):
    if c < 25:
        return 0
    if c < 37:
        return 1
    if c < 43:
        return 2

    return 3


def computeInterference(c, nConn, interMatrix):
    ret = 0.0
    for i in range(nConn):
        if c != i:
            ret += interMatrix[c][i]

    return ret


def createBigM(M_ij, nConn, interMatrix):
    for c in range(nConn):
        value = computeInterference(c, nConn, interMatrix)
        M_ij.append(value)


def convertDBMToMW(noise):
    b = noise / 10.0
    result = math.pow(10.0, b)
    return result


def loadData(
    path, receivers, senders, dataRates, SINR, interMatrix, distMatrix, M_ij,
):
    with open(path, "r") as f:
        line = f.readline()
        aux = line.split()
        nConn = int(aux[0])
        alfa = float(aux[1])
        noise = float(aux[2])
        powerSender = float(aux[3])
        nSpectrums = int(aux[4])
        specs = [0 for i in range(nSpectrums)]

        for i in range(nSpectrums):
            specs[i] = int(aux[5 + i])

        if noise != 0.0:
            noise = convertDBMToMW(noise)

        # read receivers
        # read senders
        # read data-rates
        # read SINR

        f.readline()
        for i in range(nConn):
            line = f.readline()
            aux = line.split()
            receivers[i] = [float(aux[0]), float(aux[1])]

        f.readline()
        for i in range(nConn):
            line = f.readline()
            aux = line.split()
            senders[i] = [float(aux[0]), float(aux[1])]

        f.readline()
        for i in range(nConn):
            line = f.readline()
            # gamma[i] = float(line)

        f.readline()
        for i in range(12):
            line = f.readline()
            aux = line.split()

            for j in range(4):
                dataRates[i][j] = float(aux[j])

        f.readline()
        for i in range(12):
            line = f.readline()
            aux = line.split()

            for j in range(4):
                SINR[i][j] = float(aux[j])

        convertTableToMW(SINR, SINR)

        distanceAndInterference(
            senders, receivers, interMatrix, distMatrix, powerSender, nConn, alfa,
        )

        createBigM(M_ij, nConn, interMatrix)

    return noise, powerSender, alfa, nConn


def defineModel(N, NC, MCS, DM, IM, NS, PS):
    m = gp.Model("vrbsp")

    def variables():
        x = m.addVars(nConn, nChan, 12, vtype=GRB.BINARY, name="x")
        z = m.addVars(nConn, nChan, lb=0.0, vtype=GRB.INTEGER, name="z")
        I = m.addVars(nConn, lb=0.0, vtype=GRB.CONTINUOUS, name="I")
        I_c = m.addVars(nConn, nChan, lb=0.0, vtype=GRB.CONTINUOUS, name="I_c")

        return x, z, I, I_c

    def constraints(x, z, I, I_c):
        # Constraint 1
        m.addConstrs(x.sum(i) <= 1.0 for i in range(N))

        # Constraint 2
        m.addConstrs(
            gp.quicksum(
                x[i, c2, m]
                for c2 in range(NC)
                if overlap[c1][c2] == 1
                for m in range(MCS)
            )
            == z[i, c1]
            for i in range(N)
            for c1 in range(NC)
        )

        # Constraint 3
        m.addConstrs(
            gp.quicksum(interMatrix[i][u] * z[u, c] for u in range(N) if u != i)
            == I_c[i, c]
            for i in range(N)
            for c in range(NC)
        )

        # Constraint 4
        m.addConstrs(
            I[i] >= I_c[i, c] - M_ij[i] * (1 - x[i, c, m])
            for i in range(N)
            for c in range(NC)
            for m in range(MCS)
        )

        # Constraint 5
        m.addConstrs(
            I[i] <= I_c[i, c] + M_ij[i] * (1 - x[i, c, m])
            for i in range(N)
            for c in range(NC)
            for m in range(MCS)
        )

        # Constraint 6
        m.addConstrs(
            I[i]
            <= gp.quicksum(
                x[i, c, m] * (((PS / math.pow(DM[i][i], alfa)) / SINR[m][cToB(c)]) - NS)
                for m in range(MCS)
                for c in range(NC)
            )
            for i in range(N)
        )

    def objective(x):
        m.setObjective(
            sense=GRB.MAXIMIZE,
            expr=gp.quicksum(
                dataRates[m][cToB(c)] * x[i, c, m]
                for i in range(N)
                for c in range(NC)
                for m in range(MCS)
            ),
        )

    x, w, I, I_c = variables()
    constraints(x, w, I, I_c)
    objective(x)
    return m


def run(
    nConn, interMatrix, distMatrix, dataRates, SINR, powerSender, noise, to_write,
):
    global nChan, count_inst, rst_headers
    try:
        m = defineModel(nConn, nChan, 12, distMatrix, interMatrix, noise, powerSender)
        m.write("model.lp")

        # file_log = to_write + "/log-inst" + str(inst) + ".txt"
        # m.setParam("LogFile", file_log)
        # m.setParam("LogToConsole", 0)
        m.setParam("TimeLimit", 3600)
        m.setParam("IntFeasTol", 1e-7)

        m.optimize()

        # if m.Status == GRB.OPTIMAL:
        #     file_ri = to_write + "/result_information.txt"
        #     with open(file_ri, "a") as output_re:
        #         for i in range(len(rst_headers) - 1):
        #             output_re.write(str(m.getAttr(rst_headers[i])) + " ")
        #         output_re.write(str(m.getAttr(rst_headers[len(rst_headers) - 1])))
        #         output_re.write("\n")
        #
        #     file_name = to_write + "/out-formatted" + str(inst) + ".txt"
        #     # conn, canal, bw, interference
        #     with open(file_name, "a") as f:
        #         f.write(str(m.getAttr(GRB.Attr.ObjVal)) + "\n")
        #         for i in range(nConn):
        #             for c in range(nChan):
        #                 nMCS = 12
        #                 for m in range(nMCS):
        #                     if x[i, c, m].x == 1.0:
        #                         f.write(
        #                             "%d %d %d %d %.12lf\n" % (i, c, cToB(c), m, I[i].x)
        #                         )
        #
        #     m.write("solution.sol")

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))


def loadOverlap():
    with open("./overlap.txt", "r") as f:
        idx = 0
        for line in f:
            aux = line.split(",")
            arr_aux = []
            for j in range(len(aux)):
                arr_aux.append(int(aux[j]))

            overlap.append(arr_aux)
            idx += 1


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("argument error")
        sys.exit(1)

    loadOverlap()

    # with open("result_information.txt", "a") as f:
    #     for i in range(len(rst_headers) - 1):
    #         f.write(rst_headers[i] + " ")
    #     f.write(rst_headers[len(rst_headers) - 1] + "\n")

    U_n = int(sys.argv[1])

    receivers = [[0 for i in range(2)] for _ in range(U_n)]
    senders = [[0 for i in range(2)] for _ in range(U_n)]
    dataRates = [[0 for i in range(4)] for _ in range(12)]
    SINR = [[0 for i in range(4)] for _ in range(12)]
    affectance = [[0 for i in range(U_n)] for _ in range(U_n)]
    distMatrix = [[0 for i in range(U_n)] for _ in range(U_n)]
    interMatrix = [[0 for i in range(U_n)] for _ in range(U_n)]

    M_ij = []

    inst = sys.argv[2]
    p_type = "MD-VRBSP"
    path = (
        "../instances/md-vrbsp/U_"
        + str(U_n)
        + "/"
        + p_type
        + "_U_"
        + str(U_n)
        + "_"
        + str(inst)
        + ".txt"
    )

    noise, powerSender, alfa, nConn = loadData(
        path, receivers, senders, dataRates, SINR, interMatrix, distMatrix, M_ij,
    )

    run(
        nConn,
        interMatrix,
        distMatrix,
        dataRates,
        SINR,
        powerSender,
        noise,
        sys.argv[3],
    )

################################


def defineConstraints(
    md, n, SINR, PS, NS, x, z, I, I_c,
):
    global overlap, nChan

    # Constraint 1
    md.addConstrs(x.sum(i) <= 1.0 for i in range(n))

    # Constraint 2
    md.addConstrs(
        gp.quicksum(
            x[i, c2, m]
            for c2 in range(nChan)
            if overlap[c1][c2] == 1
            for m in range(12)
        )
        == z[i, c1]
        for i in range(n)
        for c1 in range(nChan)
    )

    # Constraint 3
    md.addConstrs(
        gp.quicksum(interMatrix[i][u] * z[u, c] for u in range(n) if u != i)
        == I_c[i, c]
        for i in range(n)
        for c in range(nChan)
    )

    # Constraint 4
    md.addConstrs(
        I[i] >= I_c[i, c] - M_ij[i] * (1 - x[i, c, m])
        for i in range(n)
        for c in range(nChan)
        for m in range(12)
    )

    # Constraint 5
    md.addConstrs(
        I[i] <= I_c[i, c] + M_ij[i] * (1 - x[i, c, m])
        for i in range(n)
        for c in range(nChan)
        for m in range(12)
    )

    # Constraint 6
    md.addConstrs(
        I[i]
        <= gp.quicksum(
            x[i, c, m]
            * (((PS / math.pow(distMatrix[i][i], alfa)) / SINR[m][cToB(c)]) - NS)
            for m in range(12)
            for c in range(nChan)
        )
        for i in range(n)
    )


def defineObjectiveFunction(model, n, dataRates, x):
    global nChan

    obj = gp.quicksum(
        dataRates[m][cToB(c)] * x[i, c, m]
        for i in range(n)
        for c in range(nChan)
        for m in range(12)
    )

    model.setObjective(obj, GRB.MAXIMIZE)


def defineVariables(m, nConn):
    global nChan

    x = m.addVars(nConn, nChan, 12, lb=0.0, ub=1.0, vtype=GRB.BINARY, name="x")
    z = m.addVars(nConn, nChan, lb=0.0, ub=GRB.INFINITY, vtype=GRB.INTEGER, name="z")
    I = m.addVars(nConn, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="I")
    I_c = m.addVars(
        nConn, nChan, lb=0.0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS, name="I_c"
    )

    return x, z, I, I_c
