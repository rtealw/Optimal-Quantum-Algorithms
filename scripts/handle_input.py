from wrap_adm import runSDP, runSDPIterations
import numpy as np
import random

# Generating D
def getDAll(n):
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def getDWorstOR(n):
    return ['0' * n] + ['0' * (n-i-1) + '1' + '0' * i for i in range(n)]

# Generating E
def getEOR(D):
    return ['1' if '1' in x else '0' for x in D]

def getEParity(D):
    return [str(x.count('1') % 2) for x in D]

def getERandom(D):
    return [str(random.randint(a=0, b=1)) for x in D]

def getEFirst(D):
    return [x[0] for x in D]

def getEHalfZeros(D):
    return ['0' if x.count('0') <= len(x)//2 else '1' for x in D]

def getEAlternating(D):
    n = len(D[0])
    alter0 = '01' * (n//2) if n % 2 == 0 else '01' * (n//2) + '0'
    alter1 = '10' * (n//2) if n % 2 == 0 else '10' * (n//2) + '1'
    return ['1' if x in [alter0, alter1] else '0' for x in D]

def getESorted(D):
    return ['1' if list(x) == sorted(list(x)) else '0' for x in D]

def getERandomK(D, k):
    k = k if k <= len(D) else k % len(D)
    chosen = random.sample(D, k=k)
    return ['1' if x in chosen else '0' for x in D]

def getERandomOne(D):
    return getERandomK(D=D, k=1)

def getERandomTwo(D):
    return getERandomK(D=D, k=2)

runSDP()
runSDPIterations(iterations=5, getD=getDAll, getE=getEAlternating, start=1)