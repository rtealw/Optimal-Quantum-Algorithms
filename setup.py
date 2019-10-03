import numpy as np

D = ['00','01','10'] # inputs to Boolean function f
E = ['0', '1', '1']  # corresponding outputs to f
F = []               # Cartesian product of inputs with different outputs

n = len(D[0])                     # length of bit string
dimension = len(D) * (n + 1) + 1  # dimension of matrix

# Construct F
for i in range(len(D)):
    for j in range(len(D)):
        if E[i] != E[j]:
            F.append((i,j))

b_1s = np.ones((1,len(F))) # vector of 1s
A_1s = []
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
    A_1s.append(A_1)

b_0s = np.zeros((1, len(D))) # vector of 0s
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
    A_0s.append(A_0)
    print(A_0)
    slack_starter += 1