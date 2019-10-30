import numpy as np
import pandas as pd
from adm import solveSDP
import time
import math
import scipy.linalg
from scipy import sparse
import sys
import cProfile

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def getA0s(D, n, dimension):
    constraints = []
    for i in range(len(D)):
        slack = n * len(D) + i
        constraint = {"V" : [1, -1], "I" : [slack, dimension-1], "J" : [slack, dimension-1]}
        for j in range(n):
            coord = i * n + j
            constraint["V"] += [1]
            constraint["I"] += [coord]
            constraint["J"] += [coord]
        constraints.append(constraint)
    return constraints

def getA1s(D, n, F):
    constraints = []
    for (y,z) in F:
        constraint = {"V": [], "I" : [], "J" : []}
        for i in range(n):
            if D[y][i] != D[z][i]:
                constraint["V"] += [1]
                constraint["I"] += [n*y + i]
                constraint["J"] += [n*z + i]
        constraints.append(constraint)
    return constraints

def getConstraints(D, E):
    F = []               # Cartesian product of inputs with different outputs
    n = len(D[0])                     # length of bit string
    dimension = len(D) * (n + 1) + 1  # dimension of matrix

    # Construct F
    for i in range(len(D)):
        for j in range(len(D)):
            if E[i] != E[j]:
                F.append((i,j))
    
    constraints = []
    constraints.extend(getA0s(D=D, n=n, dimension=dimension))
    constraints.extend(getA1s(D=D, n=n, F=F))

    b_0s = np.zeros((len(D), 1), dtype = np.float32) # vector of 0s
    b_1s = np.ones((len(F), 1), dtype = np.float32) # vector of 1s
    bs = np.concatenate((b_0s, b_1s), axis=0)

    C = sparse.csr_matrix(np.zeros((dimension, dimension)))
    C[- 1, - 1] = 1
    return constraints, bs, C

D = ['000','001', '110', '110', '111'] # inputs to Boolean function f
E = ['0', '1', '1', '1', '1']  # corresponding outputs to f
D = ['11', '10', '00', '01']
E = ['1', '1', '0', '1']

def getAllBitStrings(n):
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def functionOR(bitstring):
    if '1' in bitstring:
        return '1'
    return '0'

def getE(bitstrings, function):
    return [function(bitstring) for bitstring in bitstrings]

def getORAll(n):
    D = getAllBitStrings(n)
    E = getE(D, functionOR)
    return D, E

def getORWorst(n):
    D = ['0'*n]
    E = ['0']
    for i in range(n):
        starting = '0' * (n-i-1) + '1' + '0' * i
        D.append(starting)
        E.append('1')
    return D, E

def wrapSDPSolver(D, E):
    constraints, bs, C = getConstraints(D=D, E=E)
    X = solveSDP(constraints=constraints, b=bs, C=C)
    return X[-1, -1]

def calculateSDPSolverComplexity(iterations, getDandE, filename=""):
    found_vals = []
    true_vals = []
    input_size = []
    run_time = []
    for i in range(1, iterations+1):
        print("Input size: {}".format(i))
        D, E = getDandE(i)

        print("D: {}".format(D))
        print("E: {}".format(E))

        starting_time = time.time()
        opt_val = wrapSDPSolver(D, E)
        t = time.time() - starting_time

        print("Obj. Func. Value: {}".format(opt_val))
        print("Run Time: {}  \n".format(t))

        found_vals.extend([opt_val.real])
        true_vals.extend([math.sqrt(i)])
        input_size.extend([i])
        run_time.extend([t])

    resultsDF = pd.DataFrame(data = {'Empirical': found_vals, 'Analytical': true_vals, 'n': input_size, 'RunTime': run_time})

    if filename != "":
        resultsDF.to_csv(index=False, path_or_buf= "./figures/{}.csv".format(filename))

#cProfile.run("calculateSDPSolverComplexity(20, getORWorst)", sort = "time")
#cProfile.run("calculateSDPSolverComplexity(6, getORAll)", sort = "time")
calculateSDPSolverComplexity(6, getORAll)