import numpy as np
import pandas as pd
import time
from adm import solveSDP, simplifyX
import math
import scipy.linalg
from scipy import sparse
import sys
import cProfile
from termcolor import cprint 
import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

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

    C = sparse.csr_matrix(np.zeros((dimension, dimension)))
    C[- 1, - 1] = 1
    return constraints, b, C

def getL(X, tolerance):
    '''
    Parameters:
        X : matrix X
        tolerance : how close the reconstruction has to be to X
    Returns: 
        L : matrix L such that L.H * L = X (the product of the conjugate
            transpose of L and itself is X)
    '''
    vals, vecs = np.linalg.eig(X)
    L = np.zeros(X.shape, dtype = np.complex128)
    for k in range(len(vals)):
        val = vals[k]
        vec = vecs[:,k]
        ket_k = np.zeros((len(vals),1), dtype = np.complex128)
        ket_k[k,0] = 1
        scalar = np.complex128(np.sqrt(np.absolute(val)))
        L += scalar * ket_k.dot(vec.H)

    reconstructed_X = simplifyX(np.matrix(L).H.dot(L))
    assert (np.absolute(reconstructed_X - X) < tolerance).all()
    return np.matrix(L)

def checkConstraints(L, D, E, tolerance=1e-3):
    ''' 
    Parameters:
        L : matrix such that L.H * L = X
        D : input bitstrings to f
        E : output bits of f on D
        tolerance : how closely the constraints have to be satisfied
    Returns:
        (boolean) : if the entries corresponding to the disagreeing bits
                    of every pair of inputs x,y in D sums to
                    1 if f(x) != f(y) and 0 otherwise
    '''
    n = len(D[0])
    LH = L.H
    for x_index in range(len(D)):
        x = D[x_index]
        fx = E[x_index]
        for y_index in range(x_index + 1, len(D)):
            y = D[y_index]
            fy = E[y_index]
            should_be = 1 - (fx == fy)
            summation = 0
            for i in range(len(x)):
                if x[i] != y[i]:
                    vxi = LH[x_index * n + i,:]
                    vyi = L[:,y_index * n + i]
                    summation += vxi.dot(vyi)
            assert np.absolute(should_be - summation) < tolerance
    return True

def getSpanProgram(X, D, E, tolerance=1e-4):
    little_X = X[:-len(D)-1, :-len(D)-1]
    L = getL(X=little_X, tolerance=.1)
    checkConstraints(L=L, D=D, E=E)
    n = len(D[0])
    I = []
    F0_idx = []
    for i in range(len(D)):
        if E[i] == '0':
            F0_idx.append(i)
    for x_index in F0_idx:
        x = D[x_index]
        vx = []
        for i in range(n):
            not_xi = 1 - eval(x[i])
            not_xi_vec = [0,0]
            not_xi_vec[not_xi] = 1
            vxi = np.round(np.real(L.H[n*x_index + i,:]))
            not_xi_times_vxi = np.kron(not_xi_vec, vxi)
            vx += np.asarray(not_xi_times_vxi).tolist()[0]
        I.append(vx)
    t = np.ones((len(F0_idx),1))
    #getIx(I=I, x='01', num_inputs=len(D), t=t)
    checkSpanProgram(D=D, E=E, I=I, t=t)
    return I, t

def getIx(I, x, num_inputs, t):
    I = np.array(I)
    n = len(x)
    subblock_length = n * num_inputs
    Ix = np.zeros((t.shape[0], subblock_length *  n))
    for i in range(n):
        start_index = (2 * i + eval(x[i])) * subblock_length
        current_subblock = I[:,start_index:start_index+subblock_length]
        Ix[:,subblock_length*i:subblock_length*(i+1)] = current_subblock
    return Ix

def checkSpanProgram(D, E, I, t, tolerance = 1e-4):
    I = np.array(I)
    n = len(D[0])
    subblock_length = n *len(D)
    for x_index in range(len(D)):
        Ix = getIx(I=I, x=D[x_index], num_inputs=len(D), t=t)
        t = np.matrix(t)
        linear_combo, residuals, rank, s = np.linalg.lstsq(a=Ix, b=t)
        residual = sum((np.matmul(Ix, linear_combo) - t) ** 2)[0,0]
        assert (residual < tolerance) == (E[x_index] == '1')
    return True

def wrapSDPSolver(D, E):
    constraints, b, C = getConstraints(D=D, E=E)
    X, iteration = solveSDP(constraints=constraints, b=b, C=C, accuracy=1e-6)
    V, t = getSpanProgram(X, D=D, E=E)
    return X[-1, -1], iteration

def runSDP(D, E, round_to=3):
    print("n:", len(D[0]))
    print("D:", D)
    print("E:", E)

    starting_time = time.time()
    optimal_value, num_iteration = wrapSDPSolver(D=D, E=E)
    run_time = time.time() - starting_time

    print("Optimal Query Complexity:", np.round(optimal_value, round_to))
    print("Number of Iterations:", num_iteration)
    print("Run Time:", np.round(run_time, round_to), "seconds")
    print()

    return optimal_value, num_iteration, run_time

def runSDPForN(getD, getE, n_end, n_start=1, round_to=3):
    input_sizes = []
    optimal_values = []
    run_times = []
    num_iterations = []
    for n in range(n_start, n_end + 1):
        D = getD(n=n)
        E = getE(D=D)
        optimal_value, num_iteration, run_time = runSDP(D=D, E=E, round_to=round_to)
        input_sizes += [n]
        optimal_values += [optimal_value]
        run_times += [run_time]
        num_iterations += [num_iteration]

def testSDPSolver(iterations=5, accuracy = 2):
    all_passed = True
    for n in range(1, iterations + 1):
        print("Testing SDP solver on OR for n = {}...".format(n), end=" ")
        D = [np.binary_repr(i, width=n) for i in range(2**n)]
        E = ['1' if '1' in x else '0' for x in D]
        received, iteration = wrapSDPSolver(D, E)
        expected = np.sqrt(n)
        if round(received, accuracy) == round(expected, accuracy):
            cprint("passed :)", "green")
        else:
            all_passed = False
            cprint("failed :(", "red")
            cprint("Expected: {}".format(expected), "green")
            cprint("Received: {}".format(received), "red")

    if all_passed:
        cprint("Tests passed :)", "green")
    else:
        cprint("Tests failed :(", "red")

if __name__ == '__main__':
    testSDPSolver()