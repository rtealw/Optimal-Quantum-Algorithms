import numpy as np
from quantum_query_optimizer import wrapSDPSolver

def test_OR3():
    n=3
    D = [np.binary_repr(i, width=n) for i in range(2**n)]
    E = ['1' if '1' in x else '0' for x in D]
    solution = wrapSDPSolver(D=D, E=E, run_checks=True)

    print(solution)
    assert(solution['query_complexity'] - np.sqrt(n) < 0.0001)
