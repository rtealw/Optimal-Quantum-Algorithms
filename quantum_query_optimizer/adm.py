import numpy as np
from scipy import sparse

# There are several functions in this file named after notation used in the following paper:
# http://mpc.zib.de/index.php/MPC/article/viewFile/40/20
# (Alternating direction augmented Lagrangian methods
# for semidefinite programming)

def plainA(constraints, dimension):
    '''
        Parameters:
            constraints : list of dictionaries that contains constraints
                from getContraints in constraints.py
            dimension :
                dimension of X, which is equal to dimension of C
        Returns:
            A.T : a matrix equal to the plain A in Wen et al.'s paper
    '''
    A = np.zeros((dimension**2, len(constraints)), dtype=np.float32)
    for constraint_index in range(len(constraints)):
        constraint = constraints[constraint_index]
        for entry_index in range(len(constraint["V"])):
            v = constraint["V"][entry_index]
            i = constraint["I"][entry_index]
            j = constraint["J"][entry_index]
            A[j * dimension + i, constraint_index] = v
    return A.T

def scriptA(constraints, X):
    '''
        Parameters:
            constraints : list of dictionaries that contains constraints
                from getContraints in constraints.py
            X : current matrix X which is the current solution to the SDP
        Returns: a matrix equivalent to the script A in Wen et al.'s paper (eq. 2)
    '''
    traces = []

    # Calculate trace of each constraint matrix times X
    for constraint in constraints:
        trace = 0
        V = constraint["V"]
        I = constraint["I"]
        J = constraint["J"]
        for num in range(len(V)):
            trace += X[I[num], J[num]] * V[num]
        traces.append(trace)

    return np.matrix(traces, dtype=np.float32).T # Vector of traces

def scriptAStar(A, y):
    '''
        Parameters:
            A : the result of plainA
            y : current y value as defined in Wen et al.'s paper (eq. 7)
        Returns: a matrix equivalent to the script A* in Wen et al.'s paper (p. 206)
    '''
    product = A.T.dot(y)
    dimension = int(np.sqrt(len(product)))
    return product.reshape((dimension, dimension))

def nextY(S, X, C, b, mu, pinvAAt, constraints):
    '''
        Parameters:
            S : matrix S (eq. 6b)
            X : matrix X, the current solution to the SDP (eq. 6c)
            C : objective function matrix (eq. 1)
            b : contraint vector b (eq. 1)
            mu : mu value set by user which impacts step size
            pinvAAt : pseudo-inverse of A times A transpose
            contraints : list of dictionaries that contains constraints
                from getContraints in constraints.py
        Returns: a matrix equivalent to y in Wen et al.'s paper (eq. 6a)
    '''
    matrix = -1 * pinvAAt
    vector_a = mu * (scriptA(constraints=constraints, X=X) -b)
    vector_b = scriptA(constraints=constraints, X=S-C)
    return matrix.dot(vector_a + vector_b)

def decomposeV(V):
    '''
       Parameters:
           V : matrix V as described in Wen et al.'s paper
       Returns:
           sigma_plus : positive eigenvalues of V
           Q_plus : eigenvectors of V corresponding to elements of sigma_plus
    '''
    unordered_vals, unordered_vecs = np.linalg.eigh(V)

    # Order eigenvalues and corresponding eigenvectors
    ordering = (-unordered_vals).argsort() 
    sigma = np.diag(unordered_vals[ordering])
    Q = unordered_vecs[:, ordering] 
    num_non_neg = sum(unordered_vals >= 0) # Number of non-negative eigenvalues
    
    sigma_plus = sigma[:num_non_neg, :num_non_neg]

    # Q dagger
    Q_plus = Q[:, :num_non_neg] # Get columns corresponding to positive eigenvalues
    return sigma_plus, Q_plus

def nextV(C, A, mu, X, y):
    '''
        Parameters:
            C : objective function matrix (eq. 1)
            A : the result of plainA to make scriptAStar
            mu : step size parameter
            X : matrix X, the current solution to the SDP (eq. 6c)
            y : current y value as defined in Wen et al.'s paper (eq. 7)
        Returns:
            a matrix corresponding to the current value of V
    '''
    return C - scriptAStar(A=A, y=y) - mu * X

def nextS(V):
    '''
        Parameters:
            V : matrix V as described in Wen et al.'s paper
        Returns:
           a matrix corresponding to the current value of S
    '''
    sigma_plus, Q_plus = decomposeV(V)
    first_product = np.matmul(Q_plus, sigma_plus)
    return np.matmul(first_product, Q_plus.T)

def nextX(mu, S, V):
    '''
        Parameters:
            mu : step size parameter
            S : matrix S (eq. 6b)
            V : matrix V as described in Wen et al.'s paper
        Returns:
           a matrix corresponding to X, the current solution to the SDP
    '''
    return 1/mu *(S - V)

def simplifyX(X, tolerance=1e-5):
    '''
        Parameters:
            X : matrix X, the current solution to the SDP (eq. 6c)
            tolerance : acceptable tolerance
        Returns:
            X : the original matrix X, but with values within (tolerance) of zero set to zero
    '''
    initial_shape = X.shape
    idx0 = np.absolute(X - np.zeros(initial_shape, dtype = np.float32)) < tolerance
    X[idx0] = 0
    idx1 = np.absolute(X - np.ones(initial_shape, dtype = np.float32)) < tolerance
    X[idx1] = 1
    return X

def solveSDP(constraints, b, C, accuracy=1e-5, mu=1, min_iterations=68, max_iterations=421):
    '''
        Parameters:
            contraints : list of dictionaries that contains constraints
                from getContraints in constraints.py
            tolerance : acceptable tolerance
            b : contraint vector b (eq. 1)
            C : objective function matrix (eq. 1)
            accuracy : tolerance for stopping condition. When objective function value changes by
                less than accuracy in one iteration, the algorithm stops running
            min_iterations : the minimum number of iterations to run
            max_iterations : the maximum number of iterations to run
        Returns:
            X : the solution to the specified SDP
            iteration : the number of iterations required to solve the SDP
    '''
    initial_shape = C.shape

    # Intialize values
    S = sparse.lil_matrix(np.eye(initial_shape[0], dtype=np.float32))
    X = sparse.lil_matrix(np.zeros(initial_shape, dtype=np.float32))
    old_z = X[-1, -1]
    A = plainA(constraints=constraints, dimension=initial_shape[0])
    pinvAAt = sparse.lil_matrix(np.linalg.pinv(np.matmul(A, A.T)))
    A = sparse.lil_matrix(A)

    # Iteratively solve SDP
    for iteration in range(max_iterations):
        y = nextY(S=S, X=X, C=C, b=b, mu=mu, pinvAAt = pinvAAt, constraints=constraints)
        V = nextV(C=C, A=A, mu=mu, X=X, y=y)
        S = nextS(V)
        X = nextX(mu=mu, S=S, V=V)
        X = simplifyX(X=X)

        # Check if objective value is stabilizing and stopping conditions are met
        if np.absolute(X[-1, -1] - old_z) < accuracy and iteration > min_iterations:
            break
        old_z = X[-1, -1]
    return X, iteration
