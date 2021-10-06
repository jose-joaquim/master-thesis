from pyscipopt import Model, quicksum, multidict, SCIP_PARAMSETTING


def vrbsp(n, nc, mcs):
    m = Model('vrbsp')

    x = {}
    for i in range(n):
        for j in range(nc):
            for k in range(mcs):
                x = m.addVar(vtype='B')

    
