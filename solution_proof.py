#! /usr/bin/env python3

import math
import sys

overlap = []
nChannels = 45


def gammaToBeta(gamma, dataRates, SINR, bandwidth):
    m = 8 if bandwidth == 20 else 9
    last = 0
    for i in range(-1, m, -1):
        if dataRates[i][bandwidth] >= gamma:
            last = i
        else:
            break

    return SINR[last][bandwidth]


def convertDBMToMW(noise):
    b = noise / 10.0
    result = math.pow(10.0, b)
    return result


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

        for i in range(nConnections):
            affectance.append([])
            for j in range(nConnections):
                value = powerSender / math.pow(distanceMatrix[i][j], alfa)
                affectance[i].append(value)

        beta = []
        for bandwidth in range(4):
            beta.append(gammaToBeta(gamma, dataRates, SINR, bandwidth))

        return noise, powerSender, alfa, nConnections, time_slots, beta


def cToB(c):
    if c < 25:
        return 0
    if c < 37:
        return 1
    if c < 43:
        return 2

    return 3


def readAndCheck(file, nTimeSlots, nConnections, dataRates, SINR):
    global nChannels

    connections = {}
    with open(file, "r") as f:
        for i in range(nConnections):
            aux = f.readline().split()
            t, c, id_, interference, comp_beta = aux[0], aux[1], aux[2], aux[3], aux[4]
            connections[t, c].append([id_, interference, comp_beta])

    all_fine = True
    for t in range(nTimeSlots):
        for c in range(nChannels):
            links = connections[t, c]

            # TODO: check how range() works
            for i in range(len(links) - 1, 0, -1):
                interference = 0.0
                id_i = links[i][0]
                interference_i = links[i][1]
                comp_beta_i = links[i][2]

                for j in range(i):
                    interference += affectance[links[i]][links[j]]

                # Checking interference
                if interference != interference_i:
                    print("ops")
                    all_fine = False
                    break

                bandwidth = cToB(c)
                link_sinr = affectance[i][i] / (interference + noise)

                # Descobrir b_ij minimo para transmitir com >= gamma_ij
                b_ij = gammaToBeta(gamma_ij, dataRates, SINR, bandwidth)

                if link_sinr < b_ij:
                    print("ops 2")
                    all_fine = False
                    break

                if b_ij < comp_beta_i:
                    print("ops 3")
                    all_fine = False
                    break


    return all_fine


if __name__ == "__main__":
    receivers, senders, dataRates = [[]], [[]], [[]]
    SINR, spectrums = [], []
    distanceMatrix, interferenceMatrix = [[]], [[]]
    affectance = [[]]

    inst = -1
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

    noise, power_sender, alfa, nConnections, nTimeSlots, beta = read_instance(
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


# Input: a file containing for each connection
#        the assigned channel, time-slot and computed interference.
