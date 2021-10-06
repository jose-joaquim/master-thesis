#! /usr/bin/env python3

import sys
import gurobipy as gp
from gurobipy import GRB

import mainmodule as mmod
from mainmodule import overlap
from mainmodule import nChannels
from mainmodule import rst_headers


def defineModel(N, NTS, AFF, AUXNC):
    m = gp.Model("mdvrbsp quadratic variables formulation")

    def variables():
        t = m.addVars(NTS, vtype=GRB.BINARY, name="t")
        x = m.addVars(N, len(AUXNC), NTS, vtype=GRB.BINARY, name="x")
        I_ = m.addVars(N, lb=0.0, name="I_")
        z = m.addVars(N, len(AUXNC), NTS, vtype=GRB.BINARY, name="z")
        w = m.addVars(N, N, name="w")

        return t, x, I_, z, w

    def constraints(t, x, I_, z, w):
        # 1
        m.addConstrs(t[i + 1] <= t[i] for i in range(NTS - 1))

        # 2
        m.addConstrs(
            gp.quicksum(x[i, c, t_] for c in AUXNC for t_ in range(NTS)) == 1
            for i in range(N)
        )

        # 3
        m.addConstrs(
            gp.quicksum(x[i, c, t_] for c in AUXNC) <= t[t_]
            for t_ in range(NTS)
            for i in range(N)
        )

        # 4
        m.addConstrs(
            gp.quicksum(
                x[i, c2[0], t_] for c2 in AUXNC.items() if overlap[c1[1]][c2[1]] == 1
            )
            == z[i, c1[0], t_]
            for t_ in range(NTS)
            for c1 in AUXNC.items()
            for i in range(N)
        )

        # 5
        m.addConstrs(
            w[i, u] >= x[i, c, t_] + z[u, c, t_] - 1
            for c in AUXNC
            for t_ in range(NTS)
            for u in range(N)
            for i in range(N)
            if u != i
        )

        # 6
        m.addConstrs(
            gp.quicksum(AFF[u][i] * w[i, u] for u in range(N) if u != i) == I_[i]
            for i in range(N)
        )

        # 7
        m.addConstrs(
            I_[i]
            <= gp.quicksum(
                (AFF[i][i] / beta[i][mmod.cToBIdx(c[1])] - noise) * x[i, c[0], t_]
                for t_ in range(NTS)
                for c in AUXNC.items()
            )
            for i in range(N)
        )

    def objective(t):
        m.setObjective(t.sum(), GRB.MINIMIZE)

    t, x, I_, z, w = variables()
    constraints(t, x, I_, z, w)
    objective(t)
    return m, x, I_


def optimization(
    nConnections,
    nTimeSlots,
    SINR,
    power_sender,
    noise,
    beta,
    inteferenceMatrix,
    distanceMatrix,
    affectance,
    dataRates,
    inst,
    to_write,
):
    try:
        global nChannels

        auxNC = [
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
            29,
            30,
            31,
            32,
            33,
            34,
            39,
            40,
            41,
            44,
        ]

        dicCH = {i: auxNC[i] for i in range(len(auxNC))}

        m, x, I_ = defineModel(nConnections, nTimeSlots, affectance, dicCH)

        m.write("quad_imp2.lp")
        file_log = to_write + "/log-inst" + str(inst) + ".txt"
        m.setParam("LogFile", file_log)
        m.setParam("LogToConsole", 0)
        m.setParam("TimeLimit", 60)
        m.setParam("IntFeasTol", 1e-7)

        m.optimize()

        if m.status == GRB.OPTIMAL:
            file_ri = to_write + "/result_information.txt"
            with open(file_ri, "a") as output_re:
                output_re.write(str(inst) + " ")
                for i in range(len(rst_headers) - 1):
                    output_re.write(str(m.getAttr(rst_headers[i])) + " ")
                output_re.write(str(m.getAttr(rst_headers[len(rst_headers) - 1])))
                output_re.write("\n")

            m.write("solution.sol")
            file_name = to_write + "/out-formatted" + str(inst) + ".txt"
            # conn, channel, MCS, interference
            with open(file_name, "a") as f:
                f.write(str(m.objVal) + "\n")
                for i in range(nConnections):
                    for c in range(nChannels):
                        for t in range(nTimeSlots):
                            # print("%d %d %d" % (i, c, t))
                            if x[i, c, t].x == 1.0:
                                f.write("%d %d %d %.12f\n" % (i, c, t, I_[i].x))

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("argument error")
        sys.exit(1)

    mmod.load_overlap()

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
    gamma = [0 for i in range(U_n)]

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

    print(path)

    (noise, power_sender, alfa, nConnections, nTimeSlots, beta,) = mmod.read_instance(
        path,
        receivers,
        senders,
        gamma,
        dataRates,
        SINR,
        interferenceMatrix,
        distanceMatrix,
        affectance,
    )

    ans = 0
    dic = {}
    if int(sys.argv[4]) == 1:
        ans = mmod.average_spec_qtd(nConnections, gamma, dataRates)
    elif int(sys.argv[4]) == 0:
        ans, dic = mmod.read_from_file()
        print("opa " + str(ans))
    else:
        print("invalid time-slot generator")
        sys.exit(1)

    optimization(
        nConnections,
        ans,  # nTimeSlots,
        SINR,
        power_sender,
        noise,
        beta,
        interferenceMatrix,
        distanceMatrix,
        affectance,
        dataRates,
        inst,
        sys.argv[3],
    )
