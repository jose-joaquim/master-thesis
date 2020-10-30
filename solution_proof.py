#! /usr/bin/env python3

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
