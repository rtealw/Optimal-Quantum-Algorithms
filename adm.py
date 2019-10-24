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
import scipy
import numpy.linalg as LA

#As = [[],[]]
#As[0] = np.matrix([[1,2],[3,4]])
#As[1] = np.matrix([[5,6],[7,8]])

def vec(X):
    return np.matrix(X.ravel(), dtype = np.float32).T

def mat(x):
    dimension = int(np.sqrt(len(x)))
    return x.reshape((dimension, dimension))

def plainA(As):
    A = vec(As[0])
    for i in range(1,len(As)):
        A = np.hstack((A, vec(As[i])))
    return A.T

def scriptA(As, X):
    result = []
    for matA in As:
        product = scipy.linalg.blas.sgemm(1, matA.T, X)
        result.append(np.trace(product, dtype = np.float32))
#    print(np.array(result))
    return np.matrix(result, dtype = np.float32).T

def scriptAStar(A, y):
#    print(A.shape)
#    print(y.shape)
    return mat(scipy.linalg.blas.sgemm(1,A.T, y))

# test scriptAStar

# y(k+1)
# = y(S, X) = -(scriptA scriptAStar) inverse
#              (mu (scriptA(X) - b) + scriptA(S - C))

def nextY(S, X, As, A, C, b, mu):
    matrixPart = -1 * np.linalg.pinv(scipy.linalg.blas.sgemm(1,A, A.T))
    vectorPart = mu * (scriptA(As, X) + -1 * b) + scriptA(As, S - C)
    return scipy.linalg.blas.sgemm(1, matrixPart, vectorPart)

def decomposeV(V):
    eigVals, Q = np.linalg.eig(V)  #
    ordering = (-eigVals).argsort() # puts indices in the descending order
    
    #sigma = np.diag(eigVals) #creates a big matrix sigma
 
    # to make sure our notation is correct 
    # we need to ensure that we have sigma+, sigma-
    sigma_unordered = np.diag(eigVals)
    primeV = Q.dot(sigma_unordered).dot(Q.T)
#    np.testing.assert_allclose(primeV, V, atol=1)
    assert primeV.shape == V.shape
 
    sigma = np.diag(eigVals[ordering])
    Q = Q[:, ordering] 
    
    # assert that the decomposition worked and we can reproduce V

   
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

# spectral decomposition
# S(k+1)
# = V_dag(k+1) = Q_dag sum Q_dag tranpose
# where Q sum Q tranpose
# (Q_dag Q_double_dag)((sum_plus, 0), (0, sum_minus)) (Q_dag tranpose, Q_double_dag tranpose)
def nextS(V):
    sigmaPlus, sigmaMinus, Qplus, Qminus = decomposeV(V)
    stepOne = np.matmul(Qplus, sigmaPlus)
    return np.matmul(stepOne, Qplus.T)

def nextX(mu, S, V):
    return 1/mu *(S - V)

def checkConstraints(As, bs, X, tolerance):
    result = []
    for i in range(len(As)):
        value = np.trace(scipy.linalg.blas.sgemm(1,As[i].T, X,), dtype = np.float32)
        is_satisfied = value > bs[i] - tolerance and value < bs[i] + tolerance
        difference = round(float(np.absolute(value - bs[i])), 2)
        result.extend([difference])
    #print(result)

#run this script, then run setup and proceed to the code below.
def solveSDP(As, b, C, iterations):
    mu = 1
    rho = .5
    initial_shape = np.shape(As[0])[0]
    S = np.eye(initial_shape, dtype = np.float32)
    X = np.zeros((initial_shape, initial_shape), dtype = np.float32) #np.eye(np.shape(As[0])[0])
    A = plainA(As)
    for i in range(iterations):
        y = nextY(S=S, X=X, As=As, A=A, C=C, b=b, mu=mu)
        V = nextV(C=C, A=A, mu=mu, X=X, y=y)
#        print("V")
#        print(V)
        S = nextS(V)
#        print("S")
#        print(S)
        primeX = nextX(mu=mu, S=S, V=V)
        X = (1-rho)*X + rho * primeX
        checkConstraints(As=As, bs=b, X=X, tolerance=.1)
#        print("X")
#        print(X)
#        print("S times X")
#        print(np.matmul(S, X))
        #print("i: {}".format(i+1))
        #print("obj func: {}".format(np.trace(np.matmul(C.T, X))))
        #print(scriptA(As, X))
    return X