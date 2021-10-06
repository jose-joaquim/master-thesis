#! /usr/bin/env python3

import sys
import gurobipy as gp
from gurobipy import GRB

import mainmodule as mmod
from mainmodule import overlap
from mainmodule import nChannels
from mainmodule import rst_headers


def defineVariables(model, nConnections, nTimeSlots, x_var, z_var, I_var, t_var, w_var):
    for i in range(nTimeSlots):
        name = "t[" + str(i) + "]"
        t_var[i] = model.addVar(0.0, 1.0, 1.0, GRB.BINARY, name)

    for i in range(nConnections):
        for c in range(nChannels):
            for t in range(nTimeSlots):
                name = "x[" + str(i) + "][" + str(c) + "][" + str(t) + "]"
                x_var[i, c, t] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, name)

    for i in range(nConnections):
        name = "I[" + str(i) + "]"
        I_var[i] = model.addVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, name)

    for i in range(nConnections):
        for c in range(nChannels):
            for t in range(nTimeSlots):
                name = "z[" + str(i) + "][" + str(c) + "][" + str(t) + "]"
                z_var[i, c, t] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, name)

    for i in range(nConnections):
        for j in range(nConnections):
            if i != j:
                name = "w[" + str(i) + "][" + str(j) + "]"
                w_var[i, j] = model.addVar(0.0, 1.0, 0.0, GRB.CONTINUOUS, name)


def defineConstraints(
    model,
    nConnections,
    nTimeSlots,
    x_var,
    t_var,
    I_var,
    z_var,
    w_var,
    affectance,
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
        for u in range(nConnections):
            if i != u:
                for t in range(nTimeSlots):
                    for c in range(nChannels):
                        model.addConstr(
                            w_var[i, u] >= x_var[i, c, t] + z_var[u, c, t] - 1
                        )

    # Constraint six
    for i in range(nConnections):
        expr = gp.LinExpr()
        for u in range(nConnections):
            if i != u:
                expr += affectance[u][i] * w_var[i, u]

        model.addConstr(I_var[i] == expr)

    # Constraint seven
    for i in range(nConnections):
        expr = gp.LinExpr()

        for c in range(nChannels):
            for t in range(nTimeSlots):
                value = (affectance[i][i] / beta[i][mmod.cToBIdx(c)]) - noise
                expr += value * x_var[i, c, t]

        model.addConstr(expr >= I_var[i])


def defineObjectiveFunction(model, t_var, nTimeSlots):
    obj_function = gp.LinExpr()

    for i in range(nTimeSlots):
        obj_function += t_var[i]

    model.setObjective(obj_function, GRB.MINIMIZE)


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
        model = gp.Model("md-vrbsp linear [Glover and Wolsey 1974]")
        x_var, z_var, I_var, t_var, w_var = {}, {}, {}, {}, {}

        defineVariables(
            model, nConnections, nTimeSlots, x_var, z_var, I_var, t_var, w_var
        )

        defineConstraints(
            model,
            nConnections,
            nTimeSlots,
            x_var,
            t_var,
            I_var,
            z_var,
            w_var,
            affectance,
            beta,
            noise,
        )

        defineObjectiveFunction(model, t_var, nTimeSlots)

        # model.write("model.lp")
        file_log = to_write + "/log-inst" + str(inst) + ".txt"
        model.setParam("LogFile", file_log)
        model.setParam("LogToConsole", 0)
        model.setParam("TimeLimit", 3600)
        model.optimize()

        file_ri = to_write + "/result_information.txt"
        with open(file_ri, "a") as output_re:
            output_re.write(str(inst) + " ")
            for i in range(len(rst_headers) - 1):
                output_re.write(str(model.getAttr(rst_headers[i])) + " ")
            output_re.write(str(model.getAttr(rst_headers[len(rst_headers) - 1])))
            output_re.write("\n")

        file_name = to_write + "/out-formatted" + str(inst) + ".txt"
        # conn, channel, MCS, interference
        with open(file_name, "a") as f:
            f.write(str(model.getAttr(GRB.Attr.ObjVal)) + "\n")
            for i in range(nConnections):
                for c in range(nChannels):
                    for t in range(nTimeSlots):
                        if x_var[i, c, t].getAttr("x") == 1.0:
                            f.write(
                                "%d %d %d %.12f\n" % (i, c, t, I_var[i].getAttr("x"))
                            )

        # model.write("solution.sol")

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


if __name__ == "__main__":
    if len(sys.argv) != 4:
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

    ans = mmod.determineTimeSlots(
        nConnections, interferenceMatrix, affectance, noise, beta, SINR
    )
    print(ans)

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
