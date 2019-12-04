import numpy as np
import math
from termcolor import cprint

# The process of getting input vectors from the SDP solution
# is described in Andrew Childs' notes in section 23.2 Span programs
# http://www.cs.umd.edu/~amchilds/qa/qa.pdf
# (Lecture Notes on Quantum Algorithms)

def getL(X, run_checks=True, tolerance=.1):
    '''
        Parameters:
            X : matrix X
            tolerance : how close the reconstruction has to be to X
        Returns: 
            L : matrix L such that L.H * L = X (the product of the conjugate
                transpose of L and itself is X)
    '''
    vals, vecs = np.linalg.eig(X)
    L = np.zeros(X.shape, dtype = np.complex128)
    for k in range(len(vals)):
        val = vals[k]
        vec = vecs[:,k]
        ket_k = np.zeros((len(vals),1), dtype = np.complex128)
        ket_k[k,0] = 1
        scalar = np.complex128(np.sqrt(np.absolute(val)))
        L += scalar * ket_k.dot(vec.H)

    if run_checks:
        reconstructed_X = np.matrix(L).H.dot(L)
        if not (np.absolute(reconstructed_X - X) < tolerance).all():
            cprint("Warning: The reconstruction of X from L is not close enough to X.", "yellow")
    return np.matrix(L)

def checkConstraints(L, D, E, tolerance=1e-3):
    ''' 
        Parameters:
            L : matrix such that L.H * L = X
            D : input bitstrings to f
            E : output bits of f on D
            tolerance : how closely the constraints have to be satisfied
        Returns:
            (boolean) : if the entries corresponding to the disagreeing bits
                        of every pair of inputs x,y in D sums to
                        1 if f(x) != f(y) and 0 otherwise
    '''
    is_valid = True
    n = len(D[0])
    LH = L.H
    for x_index in range(len(D)):
        x = D[x_index]
        fx = E[x_index]
        for y_index in range(x_index + 1, len(D)):
            y = D[y_index]
            fy = E[y_index]
            should_be = 1 - (fx == fy)
            summation = 0
            for i in range(len(x)):
                if x[i] != y[i]:
                    vxi = LH[x_index * n + i,:]
                    vyi = L[:,y_index * n + i]
                    summation += vxi.dot(vyi)
            if  np.absolute(should_be - summation) > tolerance:
                is_valid = False
    return is_valid

def getIx(I, x, num_inputs):
    ''' 
        Parameters:
            I : matrix of input vectors to span program
            x : element of D
            num_inputs : size of D
        Returns:
            Ix : input vectors for x
    '''
    I = np.array(I)
    num_rows = len(I)
    n = len(x)
    subblock_length = n * num_inputs
    Ix = np.zeros((num_rows, subblock_length *  n))
    # Get corresponding subblock for each bit of x
    for i in range(n):
        start_index = (2 * i + eval(x[i])) * subblock_length
        current_subblock = I[:,start_index:(start_index+subblock_length)]
        Ix[:,subblock_length*i:subblock_length*(i+1)] = current_subblock
    # Compress Ix to all non-zero vectors
    Ix = Ix[:,~np.all(Ix == 0, axis=0)]
    # If all zero vectors, return one zero vector
    if Ix.size == 0:
        Ix = np.zeros((num_rows,1))
    return Ix

def checkSpanProgram(D, E, I, t, tolerance=1e-10):
    ''' 
        Parameters:
            D : input bitstrings to f
            E : output bits on f to D
            I : input vectors of span program
            t : target vector of span program
            tolerance : how small the residual must be
        Returns:
            Ix : input vectors for x
    '''
    is_valid = True
    # Check span program determines correct output for each element of D
    for x_index in range(len(D)):
        Ix = getIx(I=I, x=D[x_index], num_inputs=len(D))
        t = np.matrix(t)
        linear_combo, residuals, rank, s = np.linalg.lstsq(a=Ix, b=t)
        closest_vector = Ix.dot(linear_combo)
        residual = np.sum(np.square(closest_vector - t))
        if (residual < tolerance) == (E[x_index] == '0'):
            is_valid = False
    return is_valid

def getSpanProgram(X, D, E, run_checks=True):
    ''' 
        Parameters:
            X : solution to SDP
            D : Boolean inputs
            E : Boolean outputs
            run_checks : whether to run checks
        Returns:
            I : input vectors of span program
            t : target vector of span program
    '''
    little_X = X[:-len(D)-1, :-len(D)-1]
    L = getL(X=little_X, run_checks=run_checks, tolerance=1e-5)

    n = len(D[0])
    I = []
    F0_idx = [] # Indices of x such that f(x) = 0
    for i in range(len(D)):
        if E[i] == '0':
            F0_idx.append(i)

    for x_index in F0_idx: # Iterate through x such that f(x) = 0
        x = D[x_index]
        vx = []
        for i in range(n): # Iterate through every bit of x
            not_xi = 1 - eval(x[i]) # Negate bit

            # Encode in vector negated bit in vector
            not_xi_vec = [0,0]
            not_xi_vec[not_xi] = 1

            # Create each v_xi
            vxi = np.real(L.H[n*x_index + i,:])
            not_xi_times_vxi = np.kron(not_xi_vec, vxi)
            vx += np.asarray(not_xi_times_vxi).tolist()[0]
        I.append(vx)

    # Set small values to zero
    I = np.matrix(I)
    I.real[abs(I.real) < 1e-5] = 0.0

    t = np.ones((len(F0_idx),1)) # Target vector dimension of F0_idx
    if run_checks:
        if not checkConstraints(L=L, D=D, E=E):
            cprint("Warning: Matrix L from SDP solution X does not satisfy constraints", "yellow")
        if not checkSpanProgram(D=D, E=E, I=I, t=t):
            cprint("Warning: Span program does not return correct output.", "yellow")
    return I, t
