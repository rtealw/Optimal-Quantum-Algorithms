import numpy as np
import time
import sys
import cProfile
from termcolor import cprint 
import warnings
from adm import solveSDP
from constraints import getConstraints
from span_program import getSpanProgram, checkSpanProgram

if not sys.warnoptions:
    warnings.simplefilter("ignore")

def wrapSDPSolver(D, E):
    starting_time = time.time()
    constraints, b, C = getConstraints(D=D, E=E)
    X, num_iteration = solveSDP(constraints=constraints, b=b, C=C, accuracy=1e-6)
    print("nQueries", X[-1,-1])
    I, t = getSpanProgram(X=X,D=D,E=E) if X[-1,-1] > 0 else 0,0
    return {
        "query_complexity" : X[-1,-1],
        "matrix_solution" : X,
        "num_iteration" : num_iteration,
        "span_vectors" : I,
        "target_vector" : t,
        "run_time" : time.time() - starting_time,
        "n_bitstring" : len(D[0])
    }

def runSDP(D, E, round_to=3):
    print("n:", len(D[0]))
    print("D:", D)
    print("E:", E)

    solution = wrapSDPSolver(D=D, E=E)

    print("Optimal Query Complexity:", np.round(solution['query_complexity'], round_to))
    print("Number of Iterations:", solution['num_iteration'])
    print("Run Time:", np.round(solution['run_time'], round_to), "seconds")
    print()

    return solution

def runSDPForN(getD, getE, n_end, n_start=1, round_to=3):
    solutions = []
    for n in range(n_start, n_end + 1):
        D = getD(n=n)
        E = getE(D=D)
        solution = runSDP(D=D, E=E, round_to=round_to)
        solutions += [solution]
    return solutions

def testSDPSolver(iterations=5, accuracy = 2):
    all_passed = True
    for n in range(1, iterations + 1):
        print("Testing SDP solver on OR for n = {}...".format(n), end=" ")
        D = [np.binary_repr(i, width=n) for i in range(2**n)]
        E = ['1' if '1' in x else '0' for x in D]
        #received, iteration, I, t = wrapSDPSolver(D, E)
        solution = wrapSDPSolver(D=D, E=E)
        expected = np.sqrt(n)
        if round(solution['query_complexity'], accuracy) == round(expected, accuracy):
            cprint("SDP solution passed :)", "green")
        else:
            all_passed = False
            cprint("SDP solution failed :(", "red")
            cprint("Expected: {}".format(expected), "green")
            cprint("Received: {}".format(received), "red")

    if all_passed:
        cprint("Tests passed :)", "green")
    else:
        cprint("Tests failed :(", "red")

def testSDPSolver(iterations=5, accuracy = 2):
    all_passed = True
    for n in range(1, iterations + 1):
        print("Testing SDP solver on OR for n = {}... \n".format(n), end=" ")
        D = [np.binary_repr(i, width=n) for i in range(2**n)]
        E = ['1' if '1' in x else '0' for x in D]
        #received, iteration, I, t = wrapSDPSolver(D, E)
        solution = wrapSDPSolver(D=D, E=E)
        expected = np.sqrt(n)

        # check optimization results
        if round(solution['query_complexity'], accuracy) == round(expected, accuracy):
            cprint("SDP solution passed :)", "green")
        else:
            all_passed = False
            cprint("SDP solution failed :(", "red")
            cprint("Expected: {}".format(expected), "green")
            cprint("Received: {}".format(received), "red")

        # check span program
        if checkSpanProgram(D, E, solution["span_vectors"], solution["target_vector"], tolerance=1e-4):
            cprint("Span solution passed :)", "green")
        else:
            cprint("Span solution failed :(", "red")

        #print(solution["span_vectors"])
    if all_passed:
        cprint("\n All Tests passed :)", "green")
    else:
        cprint("\n Tests failed :(", "red")

if __name__ == '__main__':
    testSDPSolver()