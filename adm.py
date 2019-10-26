#http://mpc.zib.de/index.php/MPC/article/viewFile/40/20
# Algorithm 1: Alternating direction augmented Lagrangian method for SDP
# Set X_0 and S_0 as semidefinite matrices
# for k in 0,1,... do
#   compute y(k+1)
#   compute V(k+1) and its eigenvalue decomposition
#       and set S(k+1) = V_dag(k+1)
#   compute X(k+1) = 1/mu (S(k+1) - V(k+1))


# V(k+1)
# = V(S, X) = C - A_curly_star(y(S,X)) - mu X
#           = C - A_curly_star(y(k+1)) - mu X

import numpy as np


#As = [[],[]]
#As[0] = np.matrix([[1,2],[3,4]])
#As[1] = np.matrix([[5,6],[7,8]])

def vec(X):
    return np.matrix(X.ravel(), dtype=np.float32).T

def mat(x):
    dimension = int(np.sqrt(len(x)))
    return x.reshape((dimension, dimension))

def plainA(As):
    A = vec(As[0])
    for i in range(1,len(As)):
        A = np.hstack((A, vec(As[i])))
    return A.T

def scriptA(As, X, A0Indices, A1Indices):
    result = []
    for current_indices in A0Indices:
        current_trace = 0
        for indice in current_indices:
            current_trace += X[indice]
        current_trace -= X[-1, -1]
        result.append(current_trace)
    for current_indices in A1Indices:
        current_trace = 0
        for indice in current_indices:
            current_trace += X[indice]
        result.append(current_trace)
    return np.matrix(result, dtype=np.float32).T

def scriptAStar(A, y):
    return mat(np.matmul(A.T, y))

def nextY(S, X, As, A, C, b, mu, pinvAAt, A0Indices, A1Indices):
    matrixPart = -1 * pinvAAt
    vectorPart = mu * (scriptA(As=As, X=X, A0Indices=A0Indices, A1Indices=A1Indices) + -1 * b) + scriptA(As=As, X=S - C, A0Indices=A0Indices, A1Indices=A1Indices)
    return np.matmul(matrixPart, vectorPart)

def decomposeV(V):
    eigVals, Q = np.linalg.eig(V) #
    ordering = (-eigVals).argsort() # puts indices in the descending order
    #sigma = np.diag(eigVals) #creates a big matrix sigma
 
    # to make sure our notation is correct 
    # we need to ensure that we have sigma+, sigma-
    sigma_unordered = np.diag(eigVals)
    primeV = Q.dot(sigma_unordered).dot(Q.T)
    assert primeV.shape == V.shape
 
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

def checkConstraints(As, bs, X, tolerance):
    result = []
    for i in range(len(As)):
        value = np.trace(np.matmul(As[i].T, X), dtype=np.float32)
        is_satisfied = value > bs[i] - tolerance and value < bs[i] + tolerance
        difference = round(float(np.absolute(value - bs[i])), 2)
        result.extend([difference])
    #print(result)

#run this script, then run setup and proceed to the code below.
def solveSDP(As, b, C, A0Indices, A1Indices, iterations):
    mu = 1
    initial_shape = np.shape(As[0])[0]
    S = np.eye(initial_shape, dtype=np.float32)
    X = np.zeros((initial_shape, initial_shape), dtype=np.float32) #np.eye(np.shape(As[0])[0])
    A = plainA(As)
    pinvAAt = np.linalg.pinv(np.matmul(A, A.T))
    for i in range(iterations):
        y = nextY(
            S=S, X=X, As=As, A=A, C=C, b=b, mu=mu, pinvAAt = pinvAAt,
            A0Indices=A0Indices, A1Indices=A1Indices
        )
        V = nextV(C=C, A=A, mu=mu, X=X, y=y)
        S = nextS(V)
        X = nextX(mu=mu, S=S, V=V)
        #checkConstraints(As=As, bs=b, X=X, tolerance=.1)
    return X