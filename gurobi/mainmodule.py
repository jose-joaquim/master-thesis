import math
import sys
import copy

overlap = []
nChannels = 45
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


def cToBIdx(c):
    if c <= 24:
        return 0
    elif c <= 36:
        return 1
    elif c <= 42:
        return 2
    else:
        return 3


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
        # distanceMatrix.append([])
        # interferenceMatrix.append([])
        X_si = receivers[i][0]
        Y_si = receivers[i][1]

        for j in range(nConnections):
            X_rj = senders[j][0]
            Y_rj = senders[j][1]

            dist = distance(X_si, Y_si, X_rj, Y_rj)
            # distanceMatrix[i].append(dist)
            distanceMatrix[i][j] = dist

            if i == j:
                interferenceMatrix[i][j] = 0.0
            else:
                value = powerSender / math.pow(dist, alfa) if dist != 0.0 else 1e9
                interferenceMatrix[i][j] = value


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
    global overlap
    with open("./overlap.txt", "r") as f:
        idx = 0
        for line in f:
            aux = line.split(",")
            arr_aux = []
            for j in range(len(aux)):
                arr_aux.append(int(aux[j]))

            overlap.append(arr_aux)
            idx += 1


def gammaToBeta(gamma, dataRates, SINR, bandwidth):
    # print(gamma, bandwidth)
    m = 11

    last = -1
    for i in range(m):
        if dataRates[i][bandwidth] >= gamma:
            last = i
            break

    if last != -1:
        return SINR[last][bandwidth]
    else:  # There is no data-rate value greater or equal than gamma
        return SINR[11][3] + 1.0


def read_instance(
    path,
    receivers,
    senders,
    gamma,
    dataRates,
    SINR,
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
        # gamma = float(aux[4])
        n_spectrums = int(aux[4])
        specs = []

        for i in range(n_spectrums):
            specs.append(int(aux[5 + i]))

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
            gamma[i] = float(line)

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

        for i in range(nConnections):
            for j in range(nConnections):
                value = powerSender / math.pow(distanceMatrix[i][j], alfa)
                affectance[i][j] = value

        beta = [[0.0 for _ in range(4)] for _ in range(nConnections)]
        for j in range(nConnections):
            for i in range(4):
                beta[j][i] = gammaToBeta(gamma[j], dataRates, SINR, i)

        return noise, powerSender, alfa, nConnections, time_slots, beta


def insertLinkInto(
    insertedLinks,
    link,
    timeslot,
    channel,
    interferenceMatrix,
    affectance,
    noise,
    beta,
    SINR,
):
    # print(
    #     "aqui ("
    #     + str(timeslot)
    #     + ", "
    #     + str(channel)
    #     + ") tem  "
    #     + str(len(insertedLinks))
    #     + " conexoes"
    # )
    if not insertedLinks:
        # print("opa")
        insertedLinks.append(link)

        m = 11
        maxChannelSINR = SINR[m][cToBIdx(channel)]

        if maxChannelSINR < beta[link][cToBIdx(channel)]:
            return False
        else:
            return True

    # Computar novas interferencias
    insertedLinks.append(link)
    for i in range(len(insertedLinks)):
        linkInterference = 0.0
        for j in range(len(insertedLinks)):
            if i != j:
                linkInterference += interferenceMatrix[j][i]

        linkSINR = sys.float_info.max
        if linkInterference != 0.0:
            linkSINR = affectance[i][i] / (linkInterference + noise)

        # print(
        #     "SINR COMPUTADO FOI "
        #     + str(linkSINR)
        #     + " beta nesse canal eh "
        #     + str(beta[link][cToBIdx(channel)])
        # )
        if linkSINR < beta[link][cToBIdx(channel)]:
            return False

    return True


def determineTimeSlots(nConnections, interferenceMatrix, affectance, noise, beta, SINR):
    global nChannels

    # for i in range(nConnections):
    #     print("beta id %d:" % i)
    #     for j in range(4):
    #         print(" %.4f" % beta[i][j], end="")
    #     print("")

    nTimeSlots = 0
    links = {}
    firstChannel = [0, 25, 37, 43]
    maxSINR = [SINR[11][idx] for idx in range(4)]

    to_sort = [(affectance[idx][idx], idx) for idx in range(nConnections)]
    to_sort.sort()
    order = [x[1] for x in to_sort]
    # print(order)
    for i in range(nConnections):
        # print("encaixando conexao " + str(i))
        ok = False
        for t in range(nTimeSlots):
            if ok:
                break

            for c in range(nChannels):
                insertedLinks = []
                # to_print = (
                #     "   testando " + str(i) + " em (" + str(t) + " " + str(c) + ")"
                # )
                # print(to_print)
                if (c, t) in links:
                    insertedLinks = copy.deepcopy(links[c, t])

                ok = insertLinkInto(
                    insertedLinks,
                    order[i],
                    t,
                    c,
                    interferenceMatrix,
                    affectance,
                    noise,
                    beta,
                    SINR,
                )

                # print(
                #     "      ->>> em (" + str(t) + ", " + str(c) + ") retornou " + str(ok)
                # )

                if ok:
                    links[c, t] = copy.deepcopy(insertedLinks)
                    # print("coloquei em [" + str(t) + ", " + str(c) + "]")
                    break

        if not ok:
            nTimeSlots += 1
            # print("      CRIEI NOVO TIME-SLOT")
            # print("coloquei em [" + str(t) + ", " + str(c) + "]")
            # Insert link in the new time slot in a compatible channel
            for idx in range(len(maxSINR)):
                if maxSINR[idx] >= beta[order[i]][idx]:
                    links[firstChannel[idx], nTimeSlots] = [order[i]]

    # for key in links:
    #     print(key, "->", links[key])

    # print("links")
    # print(links[0, 0])

    return nTimeSlots + 1
