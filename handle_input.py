from wrap_adm import wrapSDPSolver
import numpy as np
import time
import random

def runSDP(D, E, place=3):
    print("n:", len(D[0]))
    print("D:", D)
    print("E:", E)

    starting_time = time.time()
    optimal_value, num_iteration = wrapSDPSolver(D=D, E=E)
    run_time = time.time() - starting_time

    print("Optimal Query Complexity:", np.round(optimal_value, place))
    print("Number of Iterations:", num_iteration)
    print("Run Time:", np.round(run_time, place), "seconds")
    print()

    return optimal_value, num_iteration, run_time

def runSDPIterations(iterations, getD, getE, filename="", start=1, place=3):
    input_sizes = []
    optimal_values = []
    run_times = []
    num_iterations = []
    for n in range(start, iterations + 1):
        D = getD(n=n)
        E = getE(D=D)
        optimal_value, num_iteration, run_time = runSDP(D=D, E=E, place=place)
        input_sizes += [n]
        optimal_values += [optimal_value]
        run_times += [run_time]
        num_iterations += [num_iteration]

def getDAll(n):
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def getDWorstOR(n):
    return ['0' * n] + ['0' * (n-i-1) + '1' + '0' * i for i in range(n)]

def getEOR(D):
    return ['1' if '1' in x else '0' for x in D]

def getEParity(D):
    return [str(x.count('1') % 2) for x in D]

def getERandom(D):
    return [str(random.randint(a=0, b=1)) for x in D]

def getEFirst(D):
    return [x[0] for x in D]

runSDPIterations(iterations=5, getD=getDAll, getE=getERandom, start=1)