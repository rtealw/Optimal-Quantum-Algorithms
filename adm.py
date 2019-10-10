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

# S(k+1)
# = V_dag(k+1) = Q_dag sum Q_dag tranpose
# where Q sum Q tranpose
# (Q_dag Q_double_dag)((sum_plus, 0), (0, sum_minus)) (Q_dag tranpose, Q_double_dag tranpose)

import numpy as np
import np.linalg.inv as inverse

As = [[],[]]
As[0] = np.matrix([[1,2],[3,4]])
As[1] = np.matrix([[5,6],[7,8]])

def vec(X):
    return X.ravel().T

def mat(x):
    dimension = int(np.sqrt(len(x)))
    return x.reshape((dimension, dimension))

# test mat
X = As[0]

def plainA(As):
    A = vec(As[0])
    for i in range(1,len(As)):
        A = np.hstack((A, vec(As[i])))
    return A.T

# test plainA
A = plainA(As)

def scriptA(As, X):
    result = []
    for matA in As:
        product = np.matmul(matA, X)
        result.append(np.trace(product))
    return np.array(result)

# test scriptA
X = np.eye(2)
scriptA(As, X)

def scriptAStar(A, y):
    return mat(np.matmul(A.T, y))

# test scriptAStar

# y(k+1)
# = y(S, X) = -(scriptA scriptAStar) inverse
#              (mu (scriptA(X) - b) + scriptA(S - C))

def nextY(S, X, As, C, b, mu):
    A = plainA(As)
    matrixPart = -1 * inverse(np.matmul(A, A.T))
    vectorPart = mu * (scriptA(As, X) - b) + scriptA(As, S - C)
    return np.matmul(matrixPart, vectorPart)

S = np.eye(np.shape(As[0]))
X = np.eye(np.shape(As[0]))

