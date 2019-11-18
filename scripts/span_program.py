import numpy as np
import math

def getL(X, tolerance):
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

    reconstructed_X = np.matrix(L).H.dot(L)
    assert (np.absolute(reconstructed_X - X) < tolerance).all()
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
            assert np.absolute(should_be - summation) < tolerance
    return True

def getIx(I, x, num_inputs, num_rows):
    I = np.array(I)
    n = len(x)
    subblock_length = n * num_inputs
    Ix = np.zeros((num_rows, subblock_length *  n))
    for i in range(n):
        start_index = (2 * i + eval(x[i])) * subblock_length
        current_subblock = I[:,start_index:(start_index+subblock_length)]
        Ix[:,subblock_length*i:subblock_length*(i+1)] = current_subblock
    Ix = Ix[:,~np.all(Ix == 0, axis=0)]
    if Ix.size == 0:
        Ix = np.zeros((num_rows,1))
    return Ix

def checkSpanProgram(D, E, I, t, tolerance = 1e-3):
    I = np.array(I)
    print("I=", I)
    for x_index in range(len(D)):
        Ix = getIx(I=I, x=D[x_index], num_inputs=len(D), num_rows=t.shape[0])

        t = np.matrix(t)
        linear_combo, residuals, rank, s = np.linalg.lstsq(a=Ix, b=t)
        closest_vector = Ix.dot(linear_combo)
        residual = np.sum(np.square(closest_vector - t))

        print("x", D[x_index])
        print("f(x)", E[x_index])
        print("Closest_vector", closest_vector)
        print("Residual", residual)
        print("Ix=", Ix)
        print("t=", t)

        # residual < tolerance means satisfied if output is 0
        if (residual < tolerance) == (E[x_index] == '0'):
            print("x", D[x_index])
            print("f(x)", E[x_index])
            print("Closest_vector", closest_vector)
            print("Residual", residual)
            print("Ix=", Ix)
            print("t=",t)
            raise ValueError("Residual above tolerance")
            return False
    return True

def getSpanProgram(X, D, E, tolerance=1e-3, run_checks=True):
    print("entering getSpanProgram")
    little_X = X[:-len(D)-1, :-len(D)-1]
    print("little_X", little_X)
    L = getL(X=little_X, tolerance=.1)

    # print()
    assert ((np.absolute(np.matmul(L.H,L) - little_X)) < tolerance).all()

    n = len(D[0])
    I = []
    F0_idx = []
    for i in range(len(D)):
        if E[i] == '0':
            F0_idx.append(i)

    zeros = L.shape[0]
    for x_index in F0_idx: # iterate through Reichardt's F_0
        x = D[x_index]
        vx = []
        for i in range(n): # iterate through a single bit string
            not_xi = 1 - eval(x[i]) # negate this single bit

            # create bra(not xi)
            not_xi_vec = [0,0]
            not_xi_vec[not_xi] = 1

            # create each v_xi
            vxi = np.real(L.H[n*x_index + i,:])

            # print("test")
            # print("vxi", vxi)
            # for j in range(len(vxi)):
            #     print(vxi[0,j])
            #     if abs(vxi[0,j]) < tolerance:
            #         vxi[0,j] = 0
            print("vxi", vxi)
            not_xi_times_vxi = np.kron(not_xi_vec, vxi)
            vx += np.asarray(not_xi_times_vxi).tolist()[0]
        I.append(vx)
    t = np.ones((len(F0_idx),1))
    if run_checks:
        checkConstraints(L=L, D=D, E=E, tolerance=tolerance)
        checkSpanProgram(D=D, E=E, I=I, t=t, tolerance=tolerance)
    return I, t