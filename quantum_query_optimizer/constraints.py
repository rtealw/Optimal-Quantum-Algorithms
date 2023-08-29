import numpy as np
from scipy import sparse

# The constraints enforce Reichardt's
# description in Theorem 6.2 in the following paper:
# https://arxiv.org/pdf/0904.2759.pdf
# (Span programs and quantum query complexity:
# The general adversary bound is nearly tight
# for every boolean function)

def getA0s(D, dimension):
    '''
        Parameters:
            D : input bitstrings to f
            dimension : dimension of matrix X
        Returns: 
            constraints : list of dictionaries that contains constraints
                        to ensure objective value of X is the maximum
                        of all diagonal entries correspoinding to bitstrings in D
    '''
    n = len(D[0])
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

def getA1s(D, F):
    '''
        Parameters:
            D : input bitstrings to f
            F : cartesian product of y and z where y,z in {0,...,len(D)-1}
                such that f(D[y]) != f(D[z])
        Returns: 
            constraints : list of dictionaries that contains constraints
                        to ensure that entries in X corresponding to bits
                        where y and z differ sum to one for each pair y,z in F
    '''
    n = len(D[0])
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
    '''
        Parameters:
            D : input bitstrings to f
            E : output bits of f on D
        Returns: 
            constraints : list of dictionaries to ensure constraints
            b : what the sum of X times each constraint should be
            C : the matrix to extract the objective function from X
    '''
    F = []
    n = len(D[0])
    dimension = len(D) * (n + 1) + 1

    for i in range(len(D)):
        for j in range(len(D)):
            if E[i] != E[j]:
                F.append((i,j))
    
    constraints = []
    constraints.extend(getA0s(D=D, dimension=dimension))
    constraints.extend(getA1s(D=D, F=F))

    b_0 = np.zeros((len(D), 1), dtype = np.float32)
    b_1 = np.ones((len(F), 1), dtype = np.float32)
    b = np.concatenate((b_0, b_1), axis=0)

    C = sparse.lil_matrix(np.zeros((dimension, dimension)))
    C[- 1, - 1] = 1
    return constraints, b, C
