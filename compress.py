import quantum_query_optimizer as qqo
import numpy as np
from scipy.linalg import orth

def compress_dimension(vectors):
    basis = orth(vectors)
    new_vectors = basis.T @ vectors
    new_vectors /= np.linalg.norm(new_vectors, axis=0)
    new_vectors = np.where(np.isnan(new_vectors), 0, new_vectors) # nan returned if original value was 0
    return new_vectors

n = 3
D = qqo.getDAll(n=n)
E = qqo.getERandom(D=D)
solutions = qqo.runSDP(D=D, E=E)

old_span = solutions['span_vectors']
old_target = solutions['target_vector']

def compress_span(I, n, num_inputs):
    span_vectors = {}
    subblock_length = n * num_inputs
    for i in range(n):
        for xi in [0,1]:
            start_index = (2 * i + xi) * subblock_length
            current_subblock = I[:,start_index:(start_index+subblock_length)]
            span_vectors[(i, xi)] = compress_dimension(current_subblock.T).T
    return span_vectors

def check_compressed(D, E, span_vectors, t, tolerance=1e-5):
    is_valid = True
    for x_index in range(len(D)):
        x = D[x_index]
        Ix = span_vectors[(0,eval(x[0]))]
        for i in range(1,n):
            Ix = np.hstack((Ix, span_vectors[(i,eval(x[i]))]))
            

        linear_combo, residuals, rank, s = np.linalg.lstsq(a=Ix, b=t)
        closest_vector = Ix.dot(linear_combo)
        residual = np.sum(np.square(closest_vector - t))
        if (residual < tolerance) == (E[x_index] == '0'):
            print(x_index)
            is_valid = False
    return is_valid

span_vectors = compress_span(I=old_span, n=n, num_inputs=len(D))
works = check_compressed(D=D,E=E,span_vectors=span_vectors,t=old_target)
print(works)

print(old_span.shape)
numcols = sum([val.shape[1] for val in span_vectors.values()])
print(numcols)



