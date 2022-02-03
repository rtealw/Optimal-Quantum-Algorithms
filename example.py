import quantum_query_optimizer as qqo
import numpy as np
from scipy.linalg import orth

# Example 1
D = ['00', '01', '10', '11']
E = ['0', '1', '1', '1']
qqo.runSDP(D=D, E=E)

# Example 2
solutions = qqo.runSDPForN(getD=qqo.getDAll, getE=qqo.getEOR, n_end=4, n_start=4)

## Additional Example: Generate Figures
#all_solutions = qqo.runSDPForN(getD=qqo.getDAll, getE=qqo.getEOR, n_end=5, n_start=1)
#qqo.visualizeRuntime(all_solutions, title="Runtime of OR by Input Size for All Inputs", filename="figures/or_all_runtime.eps")
#qqo.visualizeComplexity(all_solutions, title="Complexity of OR by Input Size for All Inputs", filename="figures/or_all_complexity.eps")
#
#worst_solutions = qqo.runSDPForN(getD=qqo.getDWorstOR, getE=qqo.getEOR, n_end=14, n_start=1)
#qqo.visualizeRuntime(worst_solutions, title="Runtime of OR by Input Size for Worse Inputs", filename="figures/or_worst_runtime.eps")
#qqo.visualizeComplexity(worst_solutions, title="Complexity of OR by Input Size for Worse Inputs", filename="figures/or_worst_complexity.eps")

## Additional Example: Paper Figures
#all_solutions = qqo.runSDPForN(getD=qqo.getDAll, getE=qqo.getEOR, n_end=5, n_start=1)
#worst_solutions = qqo.runSDPForN(getD=qqo.getDWorstOR, getE=qqo.getEOR, n_end=5, n_start=1)
#qqo.visualizeComplexity(all_solutions, functions=[np.sqrt], labels=["Empirical", "Squareroot"], title="Complexity of OR by Input Size", filename="figures/or_complexity.eps")
#qqo.visualizeRuntime(all_solutions, more_runtimes=[worst_solutions['run_time']], labels=["All", "Worst-case"], title="Runtime of OR by Input Size", filename="figures/or_runtime.eps")
#
#parity_solutions = qqo.runSDPForN(getD=qqo.getDAll, getE=qqo.getEParity, n_end=5, n_start=1)
#qqo.visualizeComplexity(parity_solutions, functions=[lambda x:x], labels=['Empirical', 'Analytical'], title="Complexity of Parity by Input Size", filename='figures/parity_complexity.eps')

A = np.matrix([[1,0, 0], [0,1, 1], [0,1, 1], [0,1, 1]]).T # A are vectors for the span program
t = np.matrix([[1], [1], [1]]) # t is the target, it's worth noting that we always have a target vector of all ones
basis = orth(A) # create an orthogonal basis of columns of A
print("A", A)
print(np.linalg.matrix_rank(A))
print("basis", basis)

print(A.shape)
print(basis.T.shape)

# reveal compressed span program
print(basis.T @ A)
print(basis.T @ t)


def compress_span_program(vectors, target):
    basis = orth(vectors)
    new_vectors = basis.T @ vectors

    # just to make it look nice
    new_vectors /= np.linalg.norm(new_vectors, axis=0)
    new_vectors = np.where(np.isnan(new_vectors), 0, new_vectors) # nan returned if original value was 0

    new_target = basis.T @ target
    return new_vectors, new_target


old_span = solutions['span_vectors']
old_target = solutions['target_vector']

print(old_span)
print(old_target[0])
print("compressing")

new_span, new_target = compress_span_program(old_span[0], np.matrix(old_target[0]))
print(new_span.shape)
print(new_span)
print(new_target)