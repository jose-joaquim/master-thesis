#! /usr/bin/env python3

import sys

import gurobipy as gp
import math
from gurobipy import GRB

overlap = []
nChannels = 45
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
    senders,
    receivers,
    interferenceMatrix,
    distanceMatrix,
    powerSender,
    nConnections,
    alfa,
):
    for i in range(nConnections):
        X_si = receivers[i][0]
        Y_si = receivers[i][1]

        for j in range(nConnections):
            X_rj = senders[j][0]
            Y_rj = senders[j][1]

            dist = distance(X_si, Y_si, X_rj, Y_rj)
            distanceMatrix[i][j] = dist

            if i == j:
                interferenceMatrix[i][j] = 0.0
            else:
                value = powerSender / math.pow(dist, alfa) if dist != 0.0 else 1e9
                interferenceMatrix[i][j] = value


def cToB(c):
    if c < 25:
        return 0
    if c < 37:
        return 1
    if c < 43:
        return 2

    return 3


def computeInterference(c, nConnections, interferenceMatrix):
    ret = 0.0
    for i in range(nConnections):
        if c != i:
            ret += interferenceMatrix[c][i]

    return ret


def createBigM(M_ij, nConnections, interferenceMatrix):
    for c in range(nConnections):
        value = computeInterference(c, nConnections, interferenceMatrix)
        M_ij.append(value)


def convertDBMToMW(noise):
    b = noise / 10.0
    result = math.pow(10.0, b)
    return result


def loadData(
    path, receivers, senders, dataRates, SINR, interferenceMatrix, distanceMatrix, M_ij,
):
    with open(path, "r") as f:
        line = f.readline()
        aux = line.split()
        nConnections = int(aux[0])
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
        for i in range(nConnections):
            line = f.readline()
            aux = line.split()
            # receivers.append([float(aux[0]), float(aux[1])])
            receivers[i] = [float(aux[0]), float(aux[1])]

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            aux = line.split()
            # senders.append([float(aux[0]), float(aux[1])])
            senders[i] = [float(aux[0]), float(aux[1])]

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            # gamma[i] = float(line)

        f.readline()
        for i in range(12):
            line = f.readline()
            aux = line.split()

            # dataRates.append([])
            for j in range(4):
                dataRates[i][j] = float(aux[j])
                # dataRates.append(float(aux[j]))

        f.readline()
        for i in range(12):
            line = f.readline()
            aux = line.split()

            for j in range(4):
                SINR[i][j] = float(aux[j])
                # SINR[i].append(float(aux[j]))

        convertTableToMW(SINR, SINR)

        distanceAndInterference(
            senders,
            receivers,
            interferenceMatrix,
            distanceMatrix,
            powerSender,
            nConnections,
            alfa,
        )

        createBigM(M_ij, nConnections, interferenceMatrix)

    return noise, powerSender, alfa, nConnections


def defineVariables(model, nConnections, x, y, z, w, I, I_c):
    global nChannels

    for i in range(nConnections):
        for c in range(nChannels):
            name = "x[" + str(i) + "][" + str(c) + "]"
            x[i, c] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, name)

    for i in range(nConnections):
        for b in range(4):
            MCS = 12
            for m in range(MCS):
                name = "y[" + str(i) + "][" + str(b) + "][" + str(m) + "]"
                y[i, b, m] = model.addVar(0.0, 1.0, 1.0, GRB.BINARY, name)

    for i in range(nConnections):
        for c in range(nChannels):
            name = "z[" + str(i) + "][" + str(c) + "]"
            z[i, c] = model.addVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, name)

    for i in range(nConnections):
        name = "I[" + str(i) + "]"
        I[i] = model.addVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, name)

    for i in range(nConnections):
        for c in range(nChannels):
            name = "I_c[" + str(i) + "][" + str(c) + "]"
            I_c[i, c] = model.addVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, name)


def defineConstraints(
    model, nConnections, SINR, powerSender, noise, x, y, z, w, I, I_c,
):
    global overlap, nChannels
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

    # Constraint 1
    for i in range(nConnections):
        expr = gp.LinExpr()
        for c in range(nChannels):
            expr += x[i, c]

        model.addConstr(expr <= 1.0)  # TODO: verificar isso aqui

    # Constraint 2
    for b in range(4):
        for i in range(nConnections):
            expr1, expr2 = gp.LinExpr(), gp.LinExpr()
            nMCS = 12
            for m in range(nMCS):
                expr1 += y[i, b, m]

            for ch in C_b[b]:
                expr2 += x[i, ch]

            model.addConstr(expr1 <= expr2)

    # Constraint 3
    for i in range(nConnections):
        for c1 in range(nChannels):
            exp = gp.LinExpr()
            for c2 in range(nChannels):
                if overlap[c1][c2] == 1:
                    exp += x[i, c2]

            model.addConstr(exp == z[i, c1])

    # Constraint 4
    for i in range(nConnections):
        for c in range(nChannels):
            exp = gp.LinExpr()
            for u in range(nConnections):
                if i != u:
                    exp += interferenceMatrix[i][u] * z[u, c]
            model.addConstr(I_c[i, c] == exp)

    # Constraint 5
    for i in range(nConnections):
        for c in range(nChannels):
            lin_exp = gp.LinExpr()
            lin_exp = I_c[i, c] - M_ij[i] * (1 - x[i, c])
            model.addConstr(I[i] >= lin_exp)

    # Constraint 6
    for i in range(nConnections):
        for c in range(nChannels):
            lin_exp = gp.LinExpr()
            lin_exp = I_c[i, c] + M_ij[i] * (1 - x[i, c])
            model.addConstr(I[i] <= lin_exp)

    # Constraint 7s
    for i in range(nConnections):
        expr = gp.LinExpr()
        for b in range(4):
            nDataRates = 12
            for m in range(nDataRates):
                value = powerSender / math.pow(distanceMatrix[i][i], alfa)
                value /= SINR[m][b]
                value -= noise
                expr += value * y[i, b, m]

        model.addConstr(expr >= I[i])


def defineObjectiveFunction(model, nConnections, dataRates, y):
    global nChannels
    objFunction = gp.LinExpr()

    for i in range(nConnections):
        for b in range(4):
            nMCS = 12
            for m in range(nMCS):
                objFunction += dataRates[m][b] * y[i, b, m]

    model.setObjective(objFunction, GRB.MAXIMIZE)


def modelF1_v2(
    nConnections,
    interferenceMatrix,
    distanceMatrix,
    dataRates,
    SINR,
    powerSender,
    noise,
    to_write,
):
    global nChannels, count_inst, rst_headers
    try:
        # Create a new model
        model = gp.Model("vrbsp")
        # Variaveis: x, z, I, I_c, w
        x, y, z, w = {}, {}, {}, {}
        I, I_c = {}, {}
        defineVariables(model, nConnections, x, y, z, w, I, I_c)

        defineConstraints(
            model, nConnections, SINR, powerSender, noise, x, y, z, w, I, I_c,
        )
        defineObjectiveFunction(model, nConnections, dataRates, y)

        file_log = to_write + "/log-inst" + str(inst) + ".txt"
        model.setParam("LogFile", file_log)
        model.setParam("LogToConsole", 0)
        model.setParam("TimeLimit", 3600)
        model.setParam("IntFeasTol", 1e-7)
        model.optimize()

        file_ri = to_write + "/result_information.txt"
        with open(file_ri, "a") as output_re:
            for i in range(len(rst_headers) - 1):
                output_re.write(str(model.getAttr(rst_headers[i])) + " ")
            output_re.write(str(model.getAttr(rst_headers[len(rst_headers) - 1])))
            output_re.write("\n")

        file_name = to_write + "/out-formatted" + str(inst) + ".txt"
        # conn, channel, bw, interference
        with open(file_name, "a") as f:
            f.write(str(model.getAttr(GRB.Attr.ObjVal)) + "\n")
            for i in range(nConnections):
                for b in range(4):
                    nMCS = 12
                    for m in range(nMCS):
                        if y[i, b, m].getAttr("x") == 1.0:
                            normal = False
                            for c in range(nChannels):
                                if x[i, c].getAttr("x") == 1.0:
                                    normal = True
                                    f.write(
                                        "%d %d %d %d %.12f\n"
                                        % (i, c, b, m, I[i].getAttr("x"))
                                    )

                            if not normal:
                                print("PROBLEMA")
                                exit()

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


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
    distanceMatrix = [[0 for i in range(U_n)] for _ in range(U_n)]
    interferenceMatrix = [[0 for i in range(U_n)] for _ in range(U_n)]

    M_ij = []

    inst = sys.argv[2]
    print(inst)

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

    noise, powerSender, alfa, nConnections = loadData(
        path,
        receivers,
        senders,
        dataRates,
        SINR,
        interferenceMatrix,
        distanceMatrix,
        M_ij,
    )

    modelF1_v2(
        nConnections,
        interferenceMatrix,
        distanceMatrix,
        dataRates,
        SINR,
        powerSender,
        noise,
        sys.argv[3]
    )