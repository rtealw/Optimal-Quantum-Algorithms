#http://mpc.zib.de/index.php/MPC/article/viewFile/40/20

import numpy as np
from scipy import sparse
import scipy.sparse.linalg

def vec(X):
    return np.matrix(X.ravel(), dtype=np.float32).T

def plainA(As):
    A = vec(As[0])
    for i in range(1,len(As)):
        A = np.hstack((A, vec(As[i])))
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
    #sparse_unordered_vals, sparse_unordered_vecs = sparse.linalg.eigs(V)
    unordered_vals, unordered_vecs = np.linalg.eig(V)
    #print("Sparse eig vals: {}".format(sparse_unordered_vals))
    #print("Regular eig vals: {}".format(unordered_vals))
    #print("Sparse Q: {}".format(sparse_unordered_vecs))
    #print("Regular Q: {}".format(unordered_vecs))
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

def simplifyX(X, initial_shape, close_enough=1e-5):
    idx0 = np.absolute(X - np.zeros(initial_shape, dtype = np.float32)) < close_enough
    X[idx0] = 0
    idx1 = np.absolute(X - np.ones(initial_shape, dtype = np.float32)) < close_enough
    X[idx1] = 1
    return X

def getAs(constraints, dimensions):
    As = []
    for constraint in constraints:
        V = np.array(constraint["V"])
        I = np.array(constraint["I"])
        J = np.array(constraint["J"])
        A = sparse.coo_matrix((V, (I, J)), shape=dimensions).todense()
        As.append(A)
    return As

def solveSDP(constraints, b, C, accuracy=1e-5, mu=1, min_iterations=69, max_iterations=420):
    initial_shape = C.shape
    S = np.eye(initial_shape[0], dtype=np.float32)
    X = np.zeros(initial_shape, dtype=np.float32)
    old_z = X[-1, -1]

    As = getAs(constraints=constraints, dimensions=initial_shape)
    A = plainA(As)
    pinvAAt = sparse.csr_matrix(np.linalg.pinv(np.matmul(A, A.T)))
    A = sparse.csr_matrix(A)

    for iteration in range(max_iterations):
        y = nextY(S=S, X=X, C=C, b=b, mu=mu, pinvAAt = pinvAAt, constraints=constraints)
        V = nextV(C=C, A=A, mu=mu, X=X, y=y)
        S = nextS(V)
        X = nextX(mu=mu, S=S, V=V)
        X = simplifyX(X=X, initial_shape=initial_shape)

        # Check if objective value is stabilizing
        if np.absolute(X[-1, -1] - old_z) < accuracy and iteration > min_iterations:
            break
        old_z = X[-1, -1]
    print("Iterations: {}".format(iteration))
    return X