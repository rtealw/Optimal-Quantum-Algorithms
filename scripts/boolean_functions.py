import numpy as np
import random

# Generating D (set of Boolean inputs)
def getDAll(n):
    '''
        Parameters:
            n : size of bit string
        Returns: 
            D : list of all n-bit strings
    '''
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def getDWorstOR(n):
    '''
        Parameters:
            n : size of bit string
        Returns:
            D : list of all n-bit strings with at most one 1-bit
    '''
    return ['0' * n] + ['0' * (n-i-1) + '1' + '0' * i for i in range(n)]

# Generating E (corresponding output to function f)
def getEOR(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D, 0 if there are no 1s in
                and 1 otherwise
    '''
    return ['1' if '1' in x else '0' for x in D]

def getEParity(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D, 0 if there are an even number
                of 1s in x and 1 otherwise
    '''
    return [str(x.count('1') % 2) for x in D]

def getERandom(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D, randomly 0 or 1
    ''' 
    return [str(random.randint(a=0, b=1)) for x in D]

def getEFirst(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D,
                the first bit of x
    ''' 
    return [x[0] for x in D]

def getEHalfZeros(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D, 
                0 if there are more 1s than 0s and 1 otherwise
    ''' 
    return ['0' if x.count('0') < x.count('1') else '1' for x in D]

def getEAlternating(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D,
                1 if x has no consecutive 0s or 1s and 0 otherwise
    ''' 
    n = len(D[0])
    alter0 = '01' * (n//2) if n % 2 == 0 else '01' * (n//2) + '0'
    alter1 = '10' * (n//2) if n % 2 == 0 else '10' * (n//2) + '1'
    return ['1' if x in [alter0, alter1] else '0' for x in D]

def getESorted(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D,
                1 if there are no 1s before 0s in x and 0 otherwise
    ''' 
    return ['1' if list(x) == sorted(list(x)) else '0' for x in D]

def getERandomK(D, k):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D,
                1 if in random set of k elements randomly chosen from D
                and 0 otherwise
    '''
    k = k if k <= len(D) else k % len(D)
    chosen = random.sample(D, k=k)
    return ['1' if x in chosen else '0' for x in D]

def getERandomOne(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D,
                1 if in set of 1 element randomly chosen from D
                and 0 otherwise
    ''' 
    return getERandomK(D=D, k=1)

def getERandomTwo(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D,
                1 if in set of 2 elements randomly chosen from D
                and 0 otherwise
    ''' 
    return getERandomK(D=D, k=2)
