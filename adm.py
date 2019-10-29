#http://mpc.zib.de/index.php/MPC/article/viewFile/40/20
# Algorithm 1: Alternating direction augmented Lagrangian method for SDP

import numpy as np
from scipy import sparse

def vec(X):
    return np.matrix(X.ravel(), dtype=np.float32).T

def plainA(As):
    A = vec(As[0])
    for i in range(1,len(As)):
        A = np.hstack((A, vec(As[i])))
    return A.T

def mat(x):
    dimension = int(np.sqrt(len(x)))
    return x.reshape((dimension, dimension))

def scriptA(As, X, constraintAs, constraints):
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
    return mat(A.T.dot(y))

def nextY(S, X, As, A, C, b, mu, pinvAAt, constraints, constraintAs):
    matrixPart = -1 * pinvAAt
    vectorPart = mu * (scriptA(As=As, X=X, constraints=constraints, constraintAs=constraintAs) + -1 * b) + scriptA(As=As, X=S - C, constraints=constraints,  constraintAs=constraintAs)
    return matrixPart.dot(vectorPart)

def decomposeV(V):
    eigVals, Q = np.linalg.eig(V) #
    ordering = (-eigVals).argsort() # puts indices in the descending order
    #sigma = np.diag(eigVals) #creates a big matrix sigma
 
    # to make sure our notation is correct 
    # we need to ensure that we have sigma+, sigma-
    #sigma_unordered = np.diag(eigVals)
    #primeV = Q.dot(sigma_unordered).dot(Q.T)
    #assert primeV.shape == V.shape
 
    sigma = np.diag(eigVals[ordering])
    Q = Q[:, ordering] 
  
    nNonNeg = sum(eigVals >= 0) #number of non-negative eigenvalues
    
    sigmaPlus = sigma[:nNonNeg, :nNonNeg]
    sigmaMinus = sigma[nNonNeg:, nNonNeg:]
   
    #Q dagger
    Qplus = Q[:, :nNonNeg] #get columns corresponding to positive eigenvalues
    #Q double dagger    
    Qminus = Q[:, nNonNeg:] #get columns corresponding to positive eigenvalues

    return sigmaPlus, sigmaMinus, Qplus, Qminus

def nextV(C, A, mu, X, y):
    return C - scriptAStar(A=A, y=y) - mu * X

def nextS(V):
    sigmaPlus, sigmaMinus, Qplus, Qminus = decomposeV(V)
    stepOne = np.matmul(Qplus, sigmaPlus)
    return np.matmul(stepOne, Qplus.T)

def nextX(mu, S, V):
    return 1/mu *(S - V)

def simplifyX(X, initial_shape, close_enough=1e-5):
    idx0 = np.absolute(X - np.zeros((initial_shape, initial_shape), dtype = np.float32)) < close_enough
    X[idx0] = 0
    idx1 = np.absolute(X - np.ones((initial_shape, initial_shape), dtype = np.float32)) < close_enough
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

def solveSDP(As, b, C, constraints, constraintAs, accuracy=1e-4, mu=1, min_iterations=69, max_iterations=420):
    initial_shape = np.shape(As[0])[0]
    S = np.eye(initial_shape, dtype=np.float32)
    X = np.zeros((initial_shape, initial_shape), dtype=np.float32)
    old_z = X[-1, -1]

    As = getAs(constraints, C.shape)

    A = plainA(As)
    pinvAAt = sparse.csr_matrix(np.linalg.pinv(np.matmul(A, A.T)))
    sparseA = sparse.csr_matrix(A)

    for iteration in range(max_iterations):
        y = nextY(
            S=S, X=X, As=As, A=sparseA, C=C, b=b, mu=mu, pinvAAt = pinvAAt,
            constraints=constraints, constraintAs=constraintAs
        )
        V = nextV(C=C, A=sparseA, mu=mu, X=X, y=y)
        S = nextS(V)
        X = nextX(mu=mu, S=S, V=V)
        X = simplifyX(X=X, initial_shape=initial_shape)

        # Check if objective value is stagnating
        if np.absolute(X[-1, -1] - old_z) < accuracy and iteration > min_iterations:
            break
        old_z = X[-1, -1]
    print("Iterations: {}".format(iteration))
    return X