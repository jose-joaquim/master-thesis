#! /usr/bin/env python3

import gurobipy as gp
from gurobipy import GRB

import mainmodule as mmod
from mainmodule import overlap
from mainmodule import nChannels
from mainmodule import rst_headers


def defineVariables(model, nConnections, nTimeSlots, x_var, t_var, I_var):
    global nChannels
    for i in range(nTimeSlots):
        name = "t[" + str(i) + "]"
        t_var[i] = model.addVar(0.0, 1.0, 1.0, GRB.BINARY, name)

    for i in range(nConnections):
        for c in range(nChannels):
            for t in range(nTimeSlots):
                name = "x[" + str(i) + "]" + "[" + str(c) + "]" + "[" + str(t) + "]"
                x_var[i, c, t] = model.addVar(0.0, 1.0, 0.0, GRB.BINARY, name)

    for i in range(nConnections):
        name = "I[" + str(i) + "]"
        I_var[i] = model.addVar(0.0, GRB.INFINITY, 0.0, GRB.CONTINUOUS, name)


def defineObjectiveFunction(model, t_var, nTimeSlots):
    obj_function = gp.LinExpr()

    for i in range(nTimeSlots):
        obj_function += t_var[i]

    model.setObjective(obj_function, GRB.MINIMIZE)


def defineConstraints(
    model, nConnections, nTimeSlots, x_var, t_var, I_var, affectance, beta, noise,
):
    global nChannels

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
        expr = gp.LinExpr()
        for j in range(nConnections):
            if i == j:
                continue

            aff = affectance[j][i]
            for c1 in range(nChannels):
                for c2 in range(nChannels):
                    if overlap[c1][c2] == 0:
                        continue

                    for t in range(nTimeSlots):
                        expr += aff * (x_var[i, c1, t] * x_var[j, c2, t])

        model.addConstr(I_var[i] == expr)

    # Constraint four (reformulated)
    for i in range(nConnections):
        expr = gp.LinExpr()

        for c in range(nChannels):
            for t in range(nTimeSlots):
                value = (affectance[i][i] / beta[i][cToBIdx(c)]) - noise
                expr += value * x_var[i, c, t]

        model.addConstr(expr >= I_var[i])

    # # Constraint four
    # for i in range(nConnections):
    #     expr = gp.LinExpr()
    #     for c in range(nChannels):
    #         for t in range(nTimeSlots):
    #             expr += beta * x_var[i, c, t]
    #
    #     model.addConstr(affectance[i][i] / (I_var[i] + noise) >= expr)


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
        model = gp.Model("md-vrbsp [non-linear]")
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

        with open("result_information.txt", "a") as output_re:
            for i in range(len(rst_headers) - 1):
                output_re.write(str(model.getAttr(rst_headers[i])) + " ")
            output_re.write(str(model.getAttr(rst_headers[len(rst_headers) - 1])))
            output_re.write("\n")

        file_name = "out-formatted" + str(count_inst) + ".txt"
        # conn, channel, MCS, interference
        with open(file_name, "a") as f:
            f.write(str(model.getAttr(GRB.Attr.ObjVal)) + "\n")
            for i in range(nConnections):
                for c in range(nChannels):
                    for t in range(time_slots):
                        if x_var[i, c, t].getAttr("x") == 1.0:
                            f.write(
                                "%d %d %d %.12f\n" % (i, c, t, I_var[i].getAttr("x"))
                            )

        model.write("solution.sol")

    except gp.GurobiError as e:
        print("Error code " + str(e.errno) + ": " + str(e))

    except AttributeError:
        print("Encountered an attribute error")


if __name__ == "__main__":
    mmod.load_overlap()

    with open("result_information.txt", "a") as f:
        for i in range(len(rst_headers) - 1):
            f.write(rst_headers[i] + " ")
        f.write(rst_headers[len(rst_headers) - 1] + "\n")

    for idx in range(1):
        U_n = 16

        receivers = [[0 for i in range(2)] for _ in range(U_n)]
        senders = [[0 for i in range(2)] for _ in range(U_n)]
        dataRates = [[0 for i in range(4)] for _ in range(12)]
        SINR = [[0 for i in range(4)] for _ in range(12)]
        affectance = [[0 for i in range(U_n)] for _ in range(U_n)]
        distanceMatrix = [[0 for i in range(U_n)] for _ in range(U_n)]
        interferenceMatrix = [[0 for i in range(U_n)] for _ in range(U_n)]
        gamma = [0 for i in range(U_n)]

        inst = idx + 1
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

        noise, power_sender, alfa, nConnections, time_slots, beta = mmod.read_instance(
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

        # optimization(
        #     nConnections,
        #     time_slots,
        #     SINR,
        #     power_sender,
        #     noise,
        #     beta,
        #     interferenceMatrix,
        #     distanceMatrix,
        #     affectance,
        #     dataRates,
        # )
