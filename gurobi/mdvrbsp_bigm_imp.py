#! /usr/bin/env python3
# With gurobi syntatic sugars, with Big-M

import sys
import gurobipy as gp
from gurobipy import GRB

import mainmodule as mmod
from mainmodule import overlap
from mainmodule import nChannels
from mainmodule import rst_headers


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


# connections, time-slots, channels, big-M, affectance
def defineModel(N, NTS, NC, BM, AFF):
    global overlap
    m = gp.Model("mdvrbsp reduced channels")

    def variables():
        t = m.addVars(NTS, vtype=GRB.BINARY, name="t")
        x = m.addVars(N, NC, NTS, vtype=GRB.BINARY, name="x")
        I_ = m.addVars(N, lb=0.0, vtype=GRB.CONTINUOUS, name="I")
        Iij = m.addVars(N, NC, NTS, lb=0.0, vtype=GRB.CONTINUOUS, name="Iij")
        z = m.addVars(N, NC, NTS, vtype=GRB.BINARY, name="z")

        return t, x, I_, Iij, z

    def constraints(t, x, I_, Iij, z):
        # 1
        m.addConstrs(t[i + 1] <= t[i] for i in range(NTS - 1))

        # 2
        m.addConstrs(
            gp.quicksum(x[i, c, t_] for c in range(NC) for t_ in range(NTS)) == 1
            for i in range(N)
        )

        # 3
        m.addConstrs(
            gp.quicksum(x[i, c, t_] for c in range(NC)) <= t[t_]
            for t_ in range(NTS)
            for i in range(N)
        )

        # 4
        m.addConstrs(
            gp.quicksum(x[i, c2, t_] for c2 in range(NC) if overlap[c1][c2] == 1)
            == z[i, c1, t_]
            for t_ in range(NTS)
            for c1 in range(NC)
            for i in range(N)
        )

        # 5
        m.addConstrs(
            gp.quicksum(AFF[u][i] * z[u, c, t_] for u in range(N) if u != i)
            == Iij[i, c, t_]
            for t_ in range(NTS)
            for c in range(NC)
            for i in range(N)
        )

        # 6.1
        m.addConstrs(
            I_[i] >= Iij[i, c, t_] - BM[i] * (1 - x[i, c, t_])
            for t_ in range(NTS)
            for c in range(NC)
            for i in range(N)
        )

        # 6.2
        m.addConstrs(
            I_[i] <= Iij[i, c, t_] + BM[i] * (1 - x[i, c, t_])
            for t_ in range(NTS)
            for c in range(NC)
            for i in range(N)
        )

        # 7
        m.addConstrs(
            I_[i]
            <= gp.quicksum(
                (AFF[i][i] / beta[i][mmod.cToBIdx(c)] - noise) * x[i, c, t_]
                for t_ in range(NTS)
                for c in range(NC)
            )
            for i in range(N)
        )

    def objective(t):
        m.setObjective(t.sum(), GRB.MINIMIZE)

    t, x, I_, Iij, z = variables()
    constraints(t, x, I_, Iij, z)
    objective(t)
    return m


def computeBigM(nConnections, interferenceMatrix, distanceMatrix, AFF):
    ret = []
    for i in range(nConnections):
        value = 0.0
        for j in range(nConnections):
            if i != j:
                value += AFF[j][i]

        ret.append(value)

    return ret


def bwToIdx(key):
    if key <= 20:
        return 0
    elif key <= 40:
        return 1
    elif key <= 80:
        return 2
    elif key <= 160:
        return 3


def optimization(
    nConnections,
    nTimeSlots,
    dic,
    SINR,
    power_sender,
    noise,
    beta,
    inteferenceMatrix,
    distanceMatrix,
    AFF,
    dataRates,
    inst,
    to_write,
):
    try:
        # model = gp.Model("md-vrbsp linear [big-m]")
        # x_var, z_var, I_var, Iij_var, t_var = {}, {}, {}, {}, {}

        global nChannels
        bigM = computeBigM(nConnections, interferenceMatrix, distanceMatrix, AFF)
        model = defineModel(nConnections, nTimeSlots, nChannels, bigM, AFF)

        model.write("bigm_imp.lp")
        # model.optimize()

        """
        defineVariables(
            model, nConnections, nTimeSlots, x_var, z_var, I_var, Iij_var, t_var
        )

        defineConstraints(
            model,
            nConnections,
            nTimeSlots,
            x_var,
            t_var,
            I_var,
            Iij_var,
            bigM,
            z_var,
            AFF,
            beta,
            noise,
        )

        defineObjectiveFunction(model, t_var, nTimeSlots)

        # Set some values in x_var...
        used_channels = {}
        for key in dic:
            print(key, dic[key])
            bw = dic[key][0]
            bwIdx = bwToIdx(bw)

            predefined_channels = box_channels[int(key[1])][int(bwIdx)]

            print(predefined_channels)
            for el_1 in predefined_channels:
                ok = False if key in used_channels else True

                if not ok:
                    for el_2 in used_channels[key]:
                        if overlap[el_1][el_2] == 1:
                            ok = True
                            break

                if ok is True:
                    print("inseri %d" % (el_1))
                    if key in used_channels:
                        used_channels[key].append(el_1)
                    else:
                        used_channels[key] = [el_1]

                    for conns_id in dic[key][1:]:
                        print(
                            "ajustei valor de [%s, %d, %d]" % (conns_id, el_1, key[0])
                        )
                        x_var[int(conns_id), el_1, key[0]].start = 1.0
                        t_var[key[0]].start = 1.0
                    break
                else:
                    print("opa deu erro")
                    return

        model.write("model.lp")
        file_log = to_write + "/log-inst" + str(inst) + ".txt"
        model.setParam("LogFile", file_log)
        model.setParam("LogToConsole", 0)
        model.setParam("TimeLimit", 60)
        model.setParam("IntFeasTol", 1e-7)
        model.optimize()

        file_ri = to_write + "/result_information.txt"
        with open(file_ri, "a") as output_re:
            output_re.write(str(inst) + " ")
            for i in range(len(rst_headers) - 1):
                output_re.write(str(model.getAttr(rst_headers[i])) + " ")
            output_re.write(str(model.getAttr(rst_headers[len(rst_headers) - 1])))
            output_re.write("\n")

        model.write("solution.sol")
        file_name = to_write + "/out-formatted" + str(inst) + ".txt"
        # conn, channel, MCS, interference
        with open(file_name, "a") as f:
            f.write(str(model.getAttr(GRB.Attr.ObjVal)) + "\n")
            for i in range(nConnections):
                for c in range(nChannels):
                    for t in range(nTimeSlots):
                        # print("%d %d %d" % (i, c, t))
                        if x_var[i, c, t].getAttr("x") == 1.0:
                            f.write(
                                "%d %d %d %.12f\n" % (i, c, t, I_var[i].getAttr("x"))
                            )
        """

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


# instance, interval 1, path_results
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("argument error")
        sys.exit(1)

    mmod.load_overlap()

    U_n = int(sys.argv[1])

    receivers = [[0 for i in range(2)] for _ in range(U_n)]
    senders = [[0 for i in range(2)] for _ in range(U_n)]
    dataRates = [[0 for i in range(4)] for _ in range(12)]
    SINR = [[0 for i in range(4)] for _ in range(12)]
    AFF = [[0 for i in range(U_n)] for _ in range(U_n)]
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
        AFF,
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
        ans,
        dic,
        SINR,
        power_sender,
        noise,
        beta,
        interferenceMatrix,
        distanceMatrix,
        AFF,
        dataRates,
        inst,
        sys.argv[3],
    )

"""
def defineVariables(m, nConn, nTS):
    t = m.addVars(nTS, vtype=GRB.BINARY, name="t")
    x = m.addVars(nConn, nChannels, nTS, vtype=GRB.BINARY, name="x")
    I_ = m.addVars(nConn, lb=0.0, vtype=GRB.CONTINUOUS, name="I")
    Iij = m.addVars(nConn, nChannels, nTS, lb=0.0, vtype=GRB.CONTINUOUS, name="Iij")
    z = m.addVars(nConn, nChannels, nTS, vtype=GRB.BINARY, name="z")

    return t, x, I_, Iij, z


def defineConstraints(
    model,
    nConnections,
    nTimeSlots,
    x_var,
    t_var,
    I_var,
    Iij_var,
    bigM,
    z_var,
    AFF,
    beta,
    noise,
):
    # Constraint one
    for i in range(nTimeSlots - 1):
        model.addConstr(t_var[i + 1] <= t_var[i])

    # Constraint two
    for i in range(nConnections):
        expr = gp.LinExpr()
        for c in range(nChannels):
            for t in range(nTimeSlots):
                expr += x_var[i, c, t]

        model.addConstr(expr == 1.0)

    # Constraint three
    for i in range(nConnections):
        for t in range(nTimeSlots):
            expr = gp.LinExpr()
            for c in range(nChannels):
                expr += x_var[i, c, t]

            model.addConstr(expr <= t_var[t])

    # Constraint four
    for i in range(nConnections):
        for c1 in range(nChannels):
            for t in range(nTimeSlots):
                expr = gp.LinExpr()
                for c2 in range(nChannels):
                    if overlap[c1][c2] == 1:
                        expr += x_var[i, c2, t]

                model.addConstr(z_var[i, c1, t] == expr)

    # Constraint five
    for i in range(nConnections):
        for c in range(nChannels):
            for t in range(nTimeSlots):
                expr = gp.LinExpr()

                for u in range(nConnections):
                    if u != i:
                        expr += AFF[u][i] * z_var[u, c, t]

                model.addConstr(Iij_var[i, c, t] == expr)

    # Constraint six
    for i in range(nConnections):
        for c in range(nChannels):
            for t in range(nTimeSlots):
                model.addConstr(
                    I_var[i] >= Iij_var[i, c, t] - bigM[i] * (1 - x_var[i, c, t])
                )
                model.addConstr(
                    I_var[i] <= Iij_var[i, c, t] + bigM[i] * (1 - x_var[i, c, t])
                )

    # Constraint seven
    for i in range(nConnections):
        expr = gp.LinExpr()

        for c in range(nChannels):
            for t in range(nTimeSlots):
                value = (AFF[i][i] / beta[i][mmod.cToBIdx(c)]) - noise
                expr += value * x_var[i, c, t]

        model.addConstr(expr >= I_var[i])


def defineObjectiveFunction(model, t_var, nTimeSlots):
    obj_function = gp.LinExpr()

    for i in range(nTimeSlots):
        obj_function += t_var[i]

    model.setObjective(obj_function, GRB.MINIMIZE)

"""
