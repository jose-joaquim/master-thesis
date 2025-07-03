#! /usr/bin/env python3

import math
import sys
import copy
import secrets

overlap = []
nChannels = 45
rst_headers = [
    "ObjVal",
    "ObjBoundC",
    "MIPGap",
    "NumVars",
    "NumConstrs",
    "NumNZs",
    "IterCount",
    "NodeCount",
    "Runtime",
]
ch_intervals = [[0, 24], [25, 37], [38, 43], [44, 45]]

box_channels = [
    [[0, 1, 2, 3, 4, 5, 6, 7], [25, 26, 27, 28], [37, 38], [43]],
    [
        [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        [29, 30, 31, 32, 33, 34],
        [39, 40, 41],
        [44],
    ],
    [[20, 21, 22, 23, 24], [35, 36], [42]],
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
        X_si = senders[i][0]
        Y_si = senders[i][1]

        for j in range(nConnections):
            X_rj = receivers[j][0]
            Y_rj = receivers[j][1]

            dist = distance(X_si, Y_si, X_rj, Y_rj)
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
    m = 12

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
    AFF,
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
        aux_dr = []

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
            receivers[i] = [float(aux[0]), float(aux[1])]

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            aux = line.split()
            senders[i] = [float(aux[0]), float(aux[1])]

        f.readline()
        for i in range(nConnections):
            line = f.readline()
            gamma[i] = float(line)

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
                AFF[i][j] = value

        beta = [[0.0 for _ in range(4)] for _ in range(nConnections)]
        for j in range(nConnections):
            for i in range(4):
                beta[j][i] = gammaToBeta(gamma[j], dataRates, SINR, i)

        time_slots = average_spec_qtd(nConnections, gamma, dataRates)
        return noise, powerSender, alfa, nConnections, time_slots, beta


def average_spec_qtd(nConnections, gamma, dataRates):
    used_spec = 0

    for i in range(nConnections):
        bw = 0
        for j in range(4):
            if dataRates[11][j] >= gamma[i]:
                if j == 0:
                    bw = 20
                elif j == 1:
                    bw = 40
                elif j == 2:
                    bw = 80
                elif j == 3:
                    bw = 160

                break

        used_spec = used_spec + bw

    number_times = math.ceil(used_spec / 500)
    print("USEI %f ceil de %d" % (number_times, used_spec))

    return number_times


def initialize():
    load_overlap()

    U_n = 1024
    receivers = [[0 for i in range(2)] for _ in range(U_n)]
    senders = [[0 for i in range(2)] for _ in range(U_n)]
    DR = [[0 for i in range(4)] for _ in range(12)]
    SINR = [[0 for i in range(4)] for _ in range(12)]
    AFF = [[0 for i in range(U_n)] for _ in range(U_n)]
    DM = [[0 for i in range(U_n)] for _ in range(U_n)]
    IM = [[0 for i in range(U_n)] for _ in range(U_n)]
    gamma = [0 for i in range(U_n)]

    (NOI, PS, alfa, N, NTS, B) = read_instance(
        sys.argv[1],
        receivers,
        senders,
        gamma,
        DR,
        SINR,
        IM,
        DM,
        AFF,
    )

    # N :   nConnections,
    # NTS : time-slots,
    # SINR :SINR,
    # PS :  power_sender,
    # NOI : noise,
    # B :   beta,
    # IM :  intermatrix,
    # DM :  distancematrix,
    # AFF : affectance,
    # DR :  datarates,
    # gamma : inst,
    # version
    return N, NTS, SINR, PS, NOI, B, IM, DM, AFF, DR, gamma


def postProcess(m, x, I_, N, NTS, dic):
    from gurobipy import GRB

    if m.status == GRB.INFEASIBLE:
        m.computeIIS()
        m.write("conflict.ilp")
    else:
        file_ri = sys.argv[3] + "/result_information.txt"
        with open(file_ri, "a") as out_re:
            out_re.write(str(sys.argv[2]) + " ")
            for i in range(len(rst_headers) - 1):
                out_re.write(str(m.getAttr(rst_headers[i])) + " ")
            out_re.write(str(m.getAttr(rst_headers[len(rst_headers) - 1])))
            out_re.write("\n")

            # if m.status == GRB.OPTIMAL:
            #     file_name = sys.argv[3] + "/out-formatted" + str(sys.argv[2]) + ".txt"
            #     # conn, channel, MCS, interference
            #     with open(file_name, "a") as f:
            #         f.write(str(m.objVal) + "\n")
            #         for i in range(N):
            #             for c in dic:
            #                 for t in range(NTS):
            #                     if x[i, c, t].x == 1.0:
            #                         f.write("%d %d %d %.12f\n" % (i, c, t, I_[i].x))


def computeBigM(nConnections, interferenceMatrix, distanceMatrix, AFF):
    ret = []
    for i in range(nConnections):
        value = 0.0
        for j in range(nConnections):
            if i != j:
                value += AFF[j][i]

        ret.append(value)

    return ret


def opt(N, NTS, SINR, PS, NOI, B, IM, DM, AFF, DR, GMM, warm=False):
    import gurobipy as gp
    from os.path import isfile

    try:
        # 1st spec
        # auxNC = [0, 1, 2, 3, 4, 5, 6, 7, 25, 26, 27, 28, 37, 38, 43]

        # 3rd spec
        # auxNC = [20, 21, 22, 23, 24, 35, 36, 42]

        # 2nd spec
        # auxNC = [
        #     8,
        #     9,
        #     10,
        #     11,
        #     12,
        #     13,
        #     14,
        #     15,
        #     16,
        #     17,
        #     18,
        #     19,
        #     29,
        #     30,
        #     31,
        #     32,
        #     33,
        #     34,
        #     39,
        #     40,
        #     41,
        #     44,
        # ]

        auxNC = [i for i in range(45)]

        dicCH = {i: auxNC[i] for i in range(len(auxNC))}
        m = gp.Model("raw model")
        x, y, I_ = {}, {}, {}

        if sys.argv[2] == "vrbsp-ybm":
            import vrbsp_ybm_bigm as vrY

            print("vrbsp ybm")
            BM = computeBigM(N, IM, DM, AFF)
            m, x, y, I_ = vrY.defineModel(N, 4, 45, 12, overlap, DR, NOI, AFF, SINR, BM)

        elif sys.argv[2] == "disj":
            import disj_prog as disj

            print("disj ")
            BM = computeBigM(N, IM, DM, AFF)
            m, x, y = disj.defineModelDisj(
                N, 4, 45, 12, overlap, DR, NOI, AFF, SINR, BM
            )

        elif sys.argv[2] == "alternative1":
            import alternative1 as alternative1

            print("alternative1")
            BM = computeBigM(N, IM, DM, AFF)
            m, x, y = alternative1.defineAlternative1(
                N, 4, 45, 12, overlap, DR, NOI, AFF, SINR, BM
            )

        elif sys.argv[2] == "alternative2":
            import alternative2 as alternative2

            print("alternative2")
            BM = computeBigM(N, IM, DM, AFF)
            m, x, y = alternative2.defineAlternative2(
                N, 4, 45, 12, overlap, DR, NOI, AFF, SINR, BM
            )

        else:
            sys.exit(1)

        m.write("model.lp")
        m.optimize()
        m.write("sol.sol")
        print(f"\n\n obj : {m.objval:10.2f}")

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


if __name__ == "__main__":
    N, NTS, SINR, PS, NOI, B, IM, DM, AFF, DR, GMM = initialize()
    # nConnections, time-slots, SINR, power_sender, noise, beta, intermatrix, distancematrix, affectance, datarates, inst, version
    opt(N, NTS, SINR, PS, NOI, B, IM, DM, AFF, DR, GMM)
