import numpy as np
#from cvxopt import matrix
#from ubsdp import ubsdp
from adm import solveSDP

D = ['00','01','10', '11'] # inputs to Boolean function f
E = ['0', '1', '1', '1']  # corresponding outputs to f

def getA1s(F, dimension, D, n):
    A1s = []
    for (y,z) in F:
        # Construct A_1 that ensures entries corresponding
        # to bits where y and z strings different sum to 1
        A_1 = np.zeros((dimension, dimension))
        stringy = D[y]
        stringz = D[z]
        xcoord = n * y
        ycoord = n * z
        for i in range(n):
            if stringy[i] != stringz[i]:
                A_1[xcoord, ycoord] = 1
            xcoord += 1
            ycoord += 1
        A1s.append(A_1.T) # transposing because of scriptA
    return A1s

def getA0s(D, n, dimension):
    A_0s = []
    count = 0 # add constants rather than multiply
    slack_starter = n * len(D)
    for i in range(len(D)):
        # Construct A_0 that ensures input chunk
        # is less than or equal to z (using slack variables)
        A_0 = np.zeros((dimension, dimension))
        for j in range(n):
            A_0[count, count] = 1
            count+=1
        A_0[slack_starter, slack_starter] = 1
        A_0[dimension-1, dimension-1] = -1
        A_0s.append(A_0.T) # transposing because of scriptA
        slack_starter += 1
    return A_0s

def getConstraints(D, E):
    F = []               # Cartesian product of inputs with different outputs
    n = len(D[0])                     # length of bit string
    dimension = len(D) * (n + 1) + 1  # dimension of matrix

    # Construct F
    for i in range(len(D)):
        for j in range(len(D)):
            if E[i] != E[j]:
                F.append((i,j))
    
    A_1s = getA1s(F=F, dimension=dimension, D=D, n=n)

    b_1s = np.ones((len(F), 1)) # vector of 1s
    b_0s = np.zeros((len(D), 1)) # vector of 0s

    A_0s = getA0s(D=D, n=n, dimension=dimension)
    A_0s.extend(A_1s)
    As = A_0s
    bs = np.concatenate((b_1s, b_0s), axis=0)

    C = np.zeros((dimension, dimension))
    C[dimension - 1, dimension - 1] = 1
    return As, bs, C

As, bs, C = getConstraints(D=D, E=E)
X = solveSDP(As=As, b=bs, C=C, iterations=1000)
print(X)
for i in range(len(As)):
    print(np.trace(np.matmul(As[i].T, X)))
    print(bs[i])



#bs = matrix(bs)
#As = matrix(As)
#C = matrix(C)
#print(As.size)
#print(As)
#X = ubsdp(bs, As, C)