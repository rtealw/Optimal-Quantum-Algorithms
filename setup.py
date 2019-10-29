import numpy as np
import pandas as pd
from adm import solveSDP
import time
import math
import scipy.linalg
from scipy import sparse
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def getA1s(F, dimension, D, n):
    A1s = []
    A1Indices = []
    A1Values = []
    constraints = []
    for (y,z) in F:
        # Construct A_1 that ensures entries corresponding
        # to bits where y and z strings different sum to 1
        A_1 = np.zeros((dimension, dimension), dtype = np.float32)
        constraint = {"V": [], "I" : [], "J" : []}
#        stringy = D[y]
#        stringz = D[z]
        xcoord = n * y
        ycoord = n * z
        current_indices = []
        current_values = []
        for i in range(n):
#            xcoord, ycoord = n*y + i, n*z + i
            current_indices.append((xcoord, ycoord))
            current_values.append(1)
            if D[y][i] != D[z][i]:
                constraint["V"] += [1]
                constraint["I"] += [xcoord]
                constraint["J"] += [ycoord]
                A_1[xcoord, ycoord] = 1
            xcoord += 1
            ycoord += 1
        A1s.append(A_1.T) # transposing because of scriptA
        A1Indices.append(current_indices)
        A1Values.append(current_values)
        constraints.append(constraint)
#    print("A1Values: {}".format(A1Values))
#    print("A1 indices: {}".format(A1Indices))
#    print("constraintsA1: {}".format(constraints))
#    print(A1s)
    return A1s, A1Indices, A1Values, constraints

def getA0s(D, n, dimension):
    A0s = []
    A0Indices = []
    A0Values = []
    count = 0 # add constants rather than multiply
    slack_starter = n * len(D)
    slack = n * len(D)
    constraints = []
    for i in range(len(D)):
        # Construct A0 that ensures input chunk
        # is less than or equal to z (using slack variables)
        constraint = {"V" : [1, -1], "I" : [slack, dimension-1], "J" : [slack, dimension-1]}
        A0 = np.zeros((dimension, dimension), dtype = np.float32)
        current_indices = []
        current_values = []
        for j in range(n):
            #print("count: {}".format(count))
            constraint["V"] += [1]
            constraint["I"] += [count]
            constraint["J"] += [count]
            current_indices.append((count, count))
            current_values.append(1)
            A0[count, count] = 1
            if i * n + j != count:
                print(count)
                print(i * n + j)
            count+=1
        constraints.append(constraint)
        A0[slack_starter, slack_starter] = 1
        current_indices.append((slack_starter, slack_starter))
        current_values.append(1)
        A0[-1, -1] = -1
        A0s.append(A0.T) # transposing because of scriptA
        A0Indices.append(current_indices)
        current_indices.append((-1, -1))
        current_values.append(-1)
        A0Values.append(current_values)
        slack_starter += 1
        slack += 1
#    print("A0Values: {}".format(A0Values))
#    print("A0 indices: {}".format(A0Indices))
#    print("constraintsA0: {}".format(constraints))
#    print(A0s)
    return A0s, A0Indices, A0Values, constraints

def getConstraints(D, E):
    F = []               # Cartesian product of inputs with different outputs
    n = len(D[0])                     # length of bit string
    dimension = len(D) * (n + 1) + 1  # dimension of matrix

    # Construct F
    for i in range(len(D)):
        for j in range(len(D)):
            if E[i] != E[j]:
                F.append((i,j))
    
    A0s, A0Indices, A0Values, constraintsA0 = getA0s(D=D, n=n, dimension=dimension)
    A1s, A1Indices, A1Values, constraintsA1 = getA1s(F=F, dimension=dimension, D=D, n=n)
    As = []
    constraints = []
    constraints.extend(constraintsA0)
    constraints.extend(constraintsA1)
    As.extend(A0s)
    As.extend(A1s)
    constraintIndices = []
    constraintIndices.extend(A0Indices)
    constraintIndices.extend(A1Indices)
    constraintValues = []
    constraintValues.extend(A0Values)
    constraintValues.extend(A1Values)
    constraintAs = {"indices" : constraintIndices, "values" : constraintValues}

    b_0s = np.zeros((len(D), 1), dtype = np.float32) # vector of 0s
    b_1s = np.ones((len(F), 1), dtype = np.float32) # vector of 1s
    bs = np.concatenate((b_0s, b_1s), axis=0)

    C = sparse.csr_matrix(np.zeros((dimension, dimension)))
    C[- 1, - 1] = 1
    return As, bs, C, constraints, constraintAs

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

def meetsConstraints(As, bs, X, tolerance):
    for i in range(len(As)):
        value = np.trace(scipy.linalg.blas.sgemm(1, As[i].T, X))
        if value < bs[i] - tolerance or value > bs[i] + tolerance:
            print(value)
            print(bs[i])
            return False
    return True

def wrapSDPSolver(D, E):
    As, bs, C, constraints, constraintAs = getConstraints(D=D, E=E)
    X = solveSDP(As=As, b=bs, C=C, constraints = constraints, constraintAs=constraintAs)
    if not meetsConstraints(As=As, bs=bs, X=X, tolerance=1):
        raise "X does not meet constraints!"
    return X[-1, -1]

def calculateSDPSolverComplexity(iterations, getDandE, filename):
    foundVals = []
    trueVals = []
    inputSize = []
    runTime = []
    for i in range(1, iterations+1):
        print("Input size: {}".format(i))
        D, E = getDandE(i)

        print("D: {}".format(D))
        print("E: {}".format(E))

        starting_time = time.time()
        optVal = wrapSDPSolver(D, E)
        t = time.time() - starting_time

        print("Obj. Func. Value: {}".format(optVal))
        print("Run Time: {}  \n".format(t))

        foundVals.extend([optVal.real])
        trueVals.extend([math.sqrt(i)])
        inputSize.extend([i])
        runTime.extend([t])

    resultsDF = pd.DataFrame(data = {'Empirical': foundVals, 'Analytical': trueVals, 'n': inputSize, 'RunTime': runTime})

    #write ouput:
    resultsDF.to_csv(index=False, path_or_buf= "./figures/{}.csv".format(filename))

#calculateSDPSolverComplexity(20, getORWorst, "output_worst_or")
import cProfile
#cProfile.run("calculateSDPSolverComplexity(6, getORAll, 'michaelTest')", sort = "time")

calculateSDPSolverComplexity(6, getORAll, 'michaelTest')