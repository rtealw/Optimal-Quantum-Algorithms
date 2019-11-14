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

def getL(X, tolerance):
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

def checkL(L, D, E, tolerance=1e-3):
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
                xi = np.zeros((1,len(D)*n))
                yi = np.zeros((len(D) * n, 1))
                xi[0,x_index*n + i] = 1
                yi[y_index * n + i,0] = 1
                vxi = xi.dot(LH)
                vyi = L.dot(yi)
                if x[i] != y[i]:
                    summation += vxi.dot(vyi)
            assert np.absolute(should_be - summation) < tolerance

def getSpanProgram(X, D, E, tolerance=1e-4):
    little_X = X[:-len(D)-1, :-len(D)-1]
    L = getL(X=little_X, tolerance=.1)
    n = len(D[0])
    V = []
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
        V.append(vx)
    checkL(L=L, D=D, E=E)
    target = np.ones((len(F0_idx),1))
    #checkSpanProgram(D=D, E=E, V=V, target=target)
    return V, target

def checkSpanProgram(D, E, V, target, tolerance = 1e-4):
    V = np.array(V)
    n = len(D[0])
    vxi_length = n *len(D)
    for y_index in range(len(D)):
        y = D[y_index]
        I = []
        for i in range(n):
            start_index = (i*2 + eval(y[i])) *vxi_length
            current_subblock = V[:,start_index:start_index+vxi_length]
            I += current_subblock.tolist()[0]
        I = np.matrix(I)
        target = np.matrix(target)
        linear_combo, residuals, rank, s = np.linalg.lstsq(a=I, b=target)
        residual = sum((np.matmul(I, linear_combo) - target) ** 2)[0,0]
        assert (residual < tolerance) == (E[y_index] == '1')
    return True

def wrapSDPSolver(D, E):
    constraints, b, C = getConstraints(D=D, E=E)
    X, iteration = solveSDP(constraints=constraints, b=b, C=C, accuracy=1e-6)
    V, target = getSpanProgram(X, D=D, E=E)
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