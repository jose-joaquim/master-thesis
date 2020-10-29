#! /usr/bin/env python3

import gurobipy as gp
import math
import sys
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


def distance(a, b, c, d):
    return math.hypot((a - c), (b - d))


def convertDBMToMW(noise):
    b = noise / 10.0
    result = math.pow(10.0, b)
    return result


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
        distanceMatrix.append([])
        interferenceMatrix.append([])
        X_si = receivers[i][0]
        Y_si = receivers[i][1]

        for j in range(nConnections):
            X_rj = senders[j][0]
            Y_rj = senders[j][1]

            dist = distance(X_si, Y_si, X_rj, Y_rj)
            distanceMatrix[i].append(dist)

            if i == j:
                interferenceMatrix[i].append(0.0)
            else:
                value = powerSender / math.pow(dist, alfa) if dist != 0.0 else 1e9
                interferenceMatrix[i].append(value)


def convertTableToMW(SINR, SINR_):
    for i in range(len(SINR_)):
        for j in range(len(SINR_[i])):
            if SINR[i][j] != 0.0:
                b = SINR[i][j] / 10.0
                result = math.pow(10.0, b)

                SINR_[i][j] = result
            else:
                SINR_[i][j] = 0.0


def load_overlap():
    with open("./overlap.txt", "r") as f:
        idx = 0
        for line in f:
            aux = line.split(",")
            arr_aux = []
            for j in range(len(aux)):
                arr_aux.append(int(aux[j]))

            overlap.append(arr_aux)
            idx += 1


def read_instance(
    path,
    receivers,
    senders,
    dataRates,
    SINR,
    spectrums,
    interferenceMatrix,
    distanceMatrix,
    affectance,
):
    with open(path, "r") as f:
        aux = f.readline().split()
        time_slots = int(aux[0])
        nConnections = time_slots
        alfa = float(aux[1])
        noise = float(aux[2])
        powerSender = float(aux[3])
        beta = float(aux[4])
        n_spectrums = int(aux[5])
        specs = []

        for i in range(n_spectrums):
            specs.append(int(aux[6 + i]))

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
            receivers.append([float(aux[0]), float(aux[1])])

        del receivers[0]

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            aux = line.split()
            senders.append([float(aux[0]), float(aux[1])])

        del senders[0]

        f.readline()
        for i in range(12):
            line = f.readline()
            aux = line.split()

            dataRates.append([])
            for j in range(4):
                dataRates.append(float(aux[j]))

        f.readline()
        for i in range(12):
            line = f.readline()
            aux = line.split()

            SINR.append([])
            for j in range(4):
                SINR[i].append(float(aux[j]))

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

        # for i in range(nConnections):
        #     print(receivers[i])
        #
        # print("========================")
        #
        # for j in range(nConnections):
        #     print(senders[j])
        #
        # print("========================")
        #
        # for i in range(nConnections):
        #     for j in range(nConnections):
        #         print("%.4f " % distanceMatrix[i][j], end=" ")
        #     print()
        #
        # print("========================")
        #
        # for i in range(nConnections):
        #     for j in range(nConnections):
        #         print("%.4f " % interferenceMatrix[i][j], end=" ")
        #     print()
        #
        # print("========================")

        for i in range(nConnections):
            affectance.append([])
            for j in range(nConnections):
                value = powerSender / math.pow(distanceMatrix[i][j], alfa)
                affectance[i].append(value)

        # for i in range(nConnections):
        #     for j in range(nConnections):
        #         print("%.4f " % affectance[i][j], end=" ")
        #     print()

        return noise, powerSender, alfa, nConnections, time_slots, beta


def defineVariables(model, nConnections, time_slots, x_var, z_var, I_var):
    print("a")


def defineConstraints(
    model, nConnections, nTimeSlots, x_var, t_var, I_var, affectance, beta, noise,
):
    print("a")


def defineObjectiveFunction(model, z_var, time_slots):
    print("a")


def optimization(
    nConnections,
    time_slots,
    SINR,
    power_sender,
    noise,
    beta,
    inteferenceMatrix,
    distanceMatrix,
    affectance,
    dataRates,
):
    global nChannels
    try:
        model = gp.Model("md-vrbsp")
        x_var, z_var, I_var = {}, {}, {}

        defineVariables(
            model, nConnections, time_slots, x_var, z_var, I_var,
        )

        defineConstraints(
            model,
            nConnections,
            time_slots,
            x_var,
            z_var,
            I_var,
            affectance,
            beta,
            noise,
        )

        defineObjectiveFunction(model, z_var, time_slots)

        model.write("model.lp")
        model.setParam("TimeLimit", 3600)
        model.optimize()
        model.write("solution.sol")

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


if __name__ == "__main__":
    load_overlap()
    for idx in range(1):
        receivers, senders, dataRates = [[]], [[]], [[]]
        SINR, spectrums = [], []
        distanceMatrix, interferenceMatrix = [[]], [[]]
        affectance = [[]]

        inst = idx + 1
        U_n = 8
        p_type = "MDVRBSP"
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

        print(path)

        noise, power_sender, alfa, nConnections, time_slots, beta = read_instance(
            path,
            receivers,
            senders,
            dataRates,
            SINR,
            spectrums,
            interferenceMatrix,
            distanceMatrix,
            affectance,
        )

        optimization(
            nConnections,
            time_slots,
            SINR,
            power_sender,
            noise,
            beta,
            interferenceMatrix,
            distanceMatrix,
            affectance,
            dataRates,
        )
