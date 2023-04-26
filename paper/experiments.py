import quantum_query_optimizer as qqo
import numpy as np
import matplotlib.pyplot as plt
import random

def get_domain_all(n):
    '''
        Parameters:
            n : size of bit string
        Returns: 
            D : list of all n-bit strings
    '''
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def get_domain_some(n, num=32):
    '''
        Parameters:
            n : size of bit string
            num : number of bit strings to return
        Returns: 
            D : list of num n-bit strings
    '''
    num = min(num, 2**n)
    return [np.binary_repr(i, width=n) for i in random.sample(range(2**n), num)]

def get_output_random(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D, randomly 0 or 1
    ''' 
    contains_one = False
    while not contains_one:
        E = [str(random.randint(a=0, b=1)) for x in D]
        if '1' in E:
            contains_one = True
    return E

def get_output_OR(D):
    '''
        Parameters:
            D : set of input bitstrings
        Returns: 
            E : for each element x in D, 
                1 if x is a 1 in every position and 0 otherwise
    ''' 
    return ['1' if '1' in x else '0' for x in D]

def e_i(index, dimension):
    '''
        Parameters:
            i : index of basis vector
            dimension : dimension of vector space
        Returns:
            e_i : basis vector
    '''
    e_i = np.zeros(dimension)
    e_i[index] = 1
    return e_i

# this function is from qqo, but I don't think it's exported
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
            print("Warning: The reconstruction of X from L is not close enough to X.", "yellow")
    return np.matrix(L)

def solveSDP(D, E, n, min_iterations):
    '''
        Parameters:
            D : list of strings of length n
            E : list of strings of length 1
            n : the size of each string
            min_iterations : minimum number of iterations to run SDP
        Returns:
            A : query complexity
            L : matrix L such that L.H * L = X where X is the solution to the SDP
    '''
    failed = True
    while failed:
        constraints, b, C = qqo.getConstraints(D=D, E=E)
        X, num_iteration = qqo.solveSDP(constraints=constraints, b=b, C=C, accuracy=1, min_iterations=min_iterations-1)
        A = X[-1,-1]
        L = getL(X, run_checks=False)
        reconstructed_X = np.matrix(L).H.dot(L)
        if (np.absolute(reconstructed_X - X) < .01).all() and A >= 1:
            failed = False
        else:
            print('Retrying with different random D,E...')
            D = get_domain_some(n)
            E = get_output_random(D)
    return A, L

def buildPsiVxi(D, E, A, L, n):
    '''
        Parameters:
            D : list of strings of length n
            E : list of strings of length 1
            n : number of qubits
            A : query complexity
            L : matrix L such that L.H * L = X where X is the solution to the SDP
        Returns:
            psi_s : list of vectors such that f(x) = 1
            all_vxi_s : list of all vectors that form solution to SDP
    '''
    vectors = []
    all_vxi_s = {i: [] for i in range(n)}
    for index, y in enumerate(E):
        vector = np.zeros(n * len(L) * 2, dtype=np.complex128)
        for i in range(n): 
            vxi = L.H[n*index + i,:]
            vxi = np.array(np.squeeze(vxi))[0]
            all_vxi_s[i] += [vxi]
            i_vector = e_i(i, n)
            bit_vector = e_i(eval(D[index][i]), 2) if y == '1' else e_i(1 - eval(D[index][i]), 2)
            vector += np.kron(i_vector, np.kron(vxi, bit_vector))
        scale = np.sqrt(2*A)
        vector = np.insert(vector / scale, 0, 1) if y == '1' else np.insert(vector * -1 * scale, 0, 1)
        vectors += [vector / np.linalg.norm(vector)]
    return vectors, all_vxi_s

def computeError(D, E, vectors, all_vxi_s, n):
    '''
        Parameters:
            D : list of input strings of length n
            E : list of output strings of length 1
            all_vxi_s : list of all vectors that form solution to SDP
        Returns:
            max_error : maximum error in the solution to the SDP
    '''
    inner_errors = [0]
    vxi_errors = [0]
    for index1, y1 in enumerate(E):
        for index2, y2 in enumerate(E):
            if y1 != y2:
                # Check psi_x and phi_x
                inner_errors += [np.abs(np.dot(vectors[index1], vectors[index2]))]
                # Check vxi_s 
                total = 0
                for i in range(n):
                    if D[index1][i] != D[index2][i]:
                        total += np.dot(all_vxi_s[i][index1], all_vxi_s[i][index2])
                vxi_errors += [np.abs(total-1)]
    return max(inner_errors)

def printyPrettyVxi(D, n, all_vxi_s):
    '''
        Parameters:
            D : list of strings of length n
            n : number of qubits
            all_vxi_s : list of all vectors that form solution to SDP
    '''
    # Remove 0 dimensions
    for i in range(n):
        a = np.stack(all_vxi_s[i], axis=1)
        idx = np.argwhere(np.all(a[..., :] == 0, axis=1))
        a = np.delete(a, idx, axis=0)
        all_vxi_s[i] = a
    # Print remaining dimensions
    for index in range(len(D)):
        for i in range(n):
            print(f'x={D[index]}, i={i}')
            vxi = all_vxi_s[i][:,index]
            print('vxi= [', ', '.join([str(x) for x in vxi]), ']')

def fullRun(D, E, n, min_iterations, pretty_print=False, verbose=True):
    '''
        Parameters:
            D : list of strings of length n
            E : list of strings of length 1
            n : number of qubits
        Returns:
            psi : list of vectors such that f(x) = 1
    '''
    # Solve SDP numerically
    A, L = solveSDP(D, E, n, min_iterations) 

    # Build vxi vectors
    vectors, all_vxi_s = buildPsiVxi(D, E, A, L, n)
    psi_s = [vector for vector, y in zip(vectors, E) if y == '1']
    num_one_outputs = len(psi_s)
    
    # Compute error
    max_error = computeError(D, E, vectors, all_vxi_s, n)

    # Print vxi vectors
    if pretty_print:
        printyPrettyVxi(D, n, all_vxi_s)
    
    # Compute spectrum of psi_s
    psi_s = np.stack(psi_s, axis=1)
    u, s, vh = np.linalg.svd(psi_s, full_matrices=False)
    if verbose:
        print(f'Query Complexity: {A}')
        print(f'Max Error: {max_error}')
        print(f'Input Size: {len(D)}')
        print(f'Number of 1 Outputs: {num_one_outputs}')
    return s, A, max_error, num_one_outputs

# iterations vs error
def write_iterations_vs_errors(ns, num_repeats, iterations, filename):    
    for n in ns:
        for num_iter in iterations:
            for _ in range(num_repeats):
                D = get_domain_some(n)
                E = get_output_random(D)
                s, A, max_error, num_one_outputs = fullRun(D, E, n, num_iter, pretty_print=False)
                if A >= 1 and max_error < .5: # Check for failure
                    with open(filename, 'a') as f:
                        f.write(f'{n},{num_iter},{max_error}\n')
                else:
                    print('Failed :(')
    return filename

def read_iterations_vs_errors(filename):
    errors_by_iteration = {}
    for line in open(filename, 'r'):
        n, num_iter, max_error = line.split(',')
        n = int(n)
        num_iter = int(num_iter)
        max_error = float(max_error)
        if n not in errors_by_iteration:
            errors_by_iteration[n] = {}
        if num_iter not in errors_by_iteration[n]:
            errors_by_iteration[n][num_iter] = []
        errors_by_iteration[n][num_iter] += [max_error]
    return errors_by_iteration

def plot_iterations_vs_errors(errors_by_iteration):
    linestyles = [(5,(10,3)), '-', '--', '-.', ':']
    for n, linestyle in zip(errors_by_iteration, linestyles):
        iterations = []
        averages = []
        stds = []
        for num_iter in errors_by_iteration[n]:
            iterations += [num_iter]
            errors = errors_by_iteration[n][num_iter]
            averages += [np.mean(errors)]
            stds += [np.std(errors)]
        upper_confidence = [avg + std for avg, std in zip(averages, stds)]
        lower_confidence = [max(avg - std, 0) for avg, std in zip(averages, stds)]
        plt.plot(iterations, averages, linestyle=linestyle, label=f'n={n}')
        plt.fill_between(iterations, upper_confidence, lower_confidence, alpha=0.2)
    plt.title('Error vs Iterations of SDP Solver')
    plt.ylabel(r'Error ($\epsilon$)')
    plt.xlabel('Iterations of SDP Solver')
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig('iterationsVsError.pdf')
    plt.close()

# Error vs required bound
def write_error_vs_bound(ns, num_repeats, iterations, filename): 
    f = open(filename, 'a')
    for n in ns:
        for min_iterations in [iterations]*num_repeats:
            D = get_domain_some(n)
            E = get_output_random(D)
            s, A, max_error, num_one_outputs = fullRun(D, E, n, min_iterations, pretty_print=False, verbose=False)
            s_min = s[-1]
            error_cutoff = s_min / (2 * 10 * A * np.sqrt(num_one_outputs) )
            #for i in range(len(s)-1):
            #    if s[i] > 1e-5:
            #        error_candidate = (s[i] / A - s[i+1] ) * 1/np.sqrt(num_one_outputs)
            #    error_cutoff = max(error_cutoff, error_candidate)
            if s_min >= 10e-8 and max_error < .5: # Hitting precision limit or failure
                with open(filename, 'a') as f:
                    f.write(f'{n},{max_error},{error_cutoff}\n')
            else:
                print('Failed :(')
    f.close()
    return filename

def read_error_vs_bound(filename):
    errors, error_cutoffs = {}, {}
    for line in open(filename, 'r'):
        n, max_error, error_cutoff = line.split(',')
        n = int(n)
        max_error = float(max_error)
        error_cutoff = float(error_cutoff)
        if n not in errors:
            errors[n] = []
            error_cutoffs[n] = []
        errors[n] += [max_error]
        error_cutoffs[n] += [error_cutoff]
    return errors, error_cutoffs

def plot_error_vs_bound(errors, error_cutoffs):
    markers = ['v', '^', '<', '>', 'o']
    for n, marker in zip(errors, markers):
        plt.scatter(errors[n], error_cutoffs[n], marker=marker, label=f'n={n}')
    # make y scale log
    plt.yscale('log')
    plt.xscale('log')
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    x = np.linspace(0,max(ymax, xmax),100)
    # color region above line y = x green
    plt.fill_between(x, ymax, x, color='green', alpha=.1)
    plt.fill_between(x, x, 0, color='red', alpha=.1)
    plt.plot(x, x, label='x=y')   
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.title('Allowed vs Actual Error with 100 Iterations of SDP Solver')
    plt.ylabel(r'Allowed Error ($s_{\kappa}/(2 A \sqrt{c N_1})$)')
    plt.xlabel(r'Actual Error ($\epsilon$)')
    plt.legend()
    plt.tight_layout()
    plt.savefig('allowed_v_actual.pdf')
    plt.close()

ns = [5,10,15,20,25]
num_repeats = 20
iterations = list(range(10,200,20))

filename = 'iterations_vs_errors.csv'
#write_iterations_vs_errors(ns=ns, num_repeats=num_repeats, iterations=iterations, filename=filename)
errors_by_iteration = read_iterations_vs_errors(filename)
plot_iterations_vs_errors(errors_by_iteration)

filename = 'error_vs_bound.csv'
#write_error_vs_bound(ns=ns, num_repeats=num_repeats, iterations=100, filename=filename)
errors, error_cutoffs = read_error_vs_bound(filename)
plot_error_vs_bound(errors, error_cutoffs)
