#http://mpc.zib.de/index.php/MPC/article/viewFile/40/20

import numpy as np
from scipy import sparse
import scipy.sparse.linalg

def plainA(constraints, dimension):
    A = np.zeros((dimension**2, len(constraints)), dtype=np.float32)
    for constraint_index in range(len(constraints)):
        constraint = constraints[constraint_index]
        for entry_index in range(len(constraint["V"])):
            v = constraint["V"][entry_index]
            i = constraint["I"][entry_index]
            j = constraint["J"][entry_index]
            A[j * dimension + i, constraint_index] = v
    return A.T

def scriptA(constraints, X):
    traces = []
    for constraint in constraints:
        trace = 0
        V = constraint["V"]
        I = constraint["I"]
        J = constraint["J"]
        for num in range(len(V)):
            trace += X[I[num], J[num]] * V[num]
        traces.append(trace)
    return np.matrix(traces, dtype=np.float32).T

def scriptAStar(A, y):
    product = A.T.dot(y)
    dimension = int(np.sqrt(len(product)))
    return product.reshape((dimension, dimension))

def nextY(S, X, C, b, mu, pinvAAt, constraints):
    matrix = -1 * pinvAAt
    vector_a = mu * (scriptA(constraints=constraints, X=X) -b)
    vector_b = scriptA(constraints=constraints, X=S-C)
    return matrix.dot(vector_a + vector_b)

def decomposeV(V):
#    unordered_vals, unordered_vecs = sparse.linalg.eigs(V)
    unordered_vals, unordered_vecs = np.linalg.eigh(V)
    ordering = (-unordered_vals).argsort() 
    sigma = np.diag(unordered_vals[ordering])
    Q = unordered_vecs[:, ordering] 
    num_non_neg = sum(unordered_vals >= 0) #number of non-negative eigenvalues
    
    sigma_plus = sigma[:num_non_neg, :num_non_neg]
    #Q dagger
    Q_plus = Q[:, :num_non_neg] #get columns corresponding to positive eigenvalues
    return sigma_plus, Q_plus

def nextV(C, A, mu, X, y):
    return C - scriptAStar(A=A, y=y) - mu * X

def nextS(V):
    sigma_plus, Q_plus = decomposeV(V)
    first_product = np.matmul(Q_plus, sigma_plus)
    return np.matmul(first_product, Q_plus.T)

def nextX(mu, S, V):
    return 1/mu *(S - V)

def simplifyX(X, close_enough=1e-5):
    initial_shape = X.shape
    idx0 = np.absolute(X - np.zeros(initial_shape, dtype = np.float32)) < close_enough
    X[idx0] = 0
    idx1 = np.absolute(X - np.ones(initial_shape, dtype = np.float32)) < close_enough
    X[idx1] = 1
    return X

def solveSDP(constraints, b, C, accuracy=1e-5, mu=1, min_iterations=68, max_iterations=421):
    initial_shape = C.shape
    S = sparse.csr_matrix(np.eye(initial_shape[0], dtype=np.float32))
    X = sparse.csr_matrix(np.zeros(initial_shape, dtype=np.float32))
    old_z = X[-1, -1]

    A = plainA(constraints=constraints, dimension=initial_shape[0])
    pinvAAt = sparse.csr_matrix(np.linalg.pinv(np.matmul(A, A.T)))
    A = sparse.csr_matrix(A)

    for iteration in range(max_iterations):
        y = nextY(S=S, X=X, C=C, b=b, mu=mu, pinvAAt = pinvAAt, constraints=constraints)
        V = nextV(C=C, A=A, mu=mu, X=X, y=y)
        S = nextS(V)
        X = nextX(mu=mu, S=S, V=V)
        X = simplifyX(X=X)

        # Check if objective value is stabilizing
        if np.absolute(X[-1, -1] - old_z) < accuracy and iteration > min_iterations:
            break
        old_z = X[-1, -1]
    return X, iteration