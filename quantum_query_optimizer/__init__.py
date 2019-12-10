import numpy as np
import time
import sys
import cProfile
import warnings

import sys
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)

from adm import solveSDP
from constraints import getConstraints
from span_program import getSpanProgram

# Pass through to user of this file
from visualize import *
from boolean_functions import *

if not sys.warnoptions:
    warnings.simplefilter("ignore")

def wrapSDPSolver(D, E, run_checks=True):
    '''
        Parameters:
            D : Boolean inputs
            E : Boolean outputs
            run_checks : whether to run checks on span program
        Returns:
            solution : dictionary with metadata from call to SDP
                query_complexity : optimal quantum query complexity of function from D to E
                matrix_solution : solution X to SDP
                num_iteration : number of iterations solveSDP ran for
                span_vectors : input vectors of span program
                target_vector : target vector of span program
                run_time : time (seconds) to run whole function
                n_bitstring : length of bitstring
    '''
    starting_time = time.time()
    constraints, b, C = getConstraints(D=D, E=E)
    X, num_iteration = solveSDP(constraints=constraints, b=b, C=C, accuracy=1e-6)
    I, t = getSpanProgram(X=X,D=D,E=E,run_checks=run_checks) if X[-1,-1] > 0 else (0,0)
    return {
        "query_complexity" : X[-1,-1],
        "matrix_solution" : X,
        "num_iteration" : num_iteration,
        "span_vectors" : I,
        "target_vector" : t,
        "run_time" : time.time() - starting_time,
        "n_bitstring" : len(D[0])
    }

def runSDP(D, E, print_output=True, round_to=3, run_checks=True):
    '''
        Parameters:
            D : Boolean inputs
            E : Boolean outputs
            print_output : whether to print output of wrapSPDSolver
            round_to : accuracy of query complexity and time (seconds)
            run_checks : whether to run checks on span program
        Returns:
            solution : dictionary with metadata from call to SDP
                query_complexity : optimal quantum query complexity of function from D to E
                matrix_solution : solution X to SDP
                num_iteration : number of iterations solveSDP ran for
                span_vectors : input vectors of span program
                target_vector : target vector of span program
                run_time : time (seconds) to run whole function
                n_bitstring : length of bitstring
    '''
    if print_output:
        print("n:", len(D[0]))
        print("D:", D)
        print("E:", E)

    solution = wrapSDPSolver(D=D, E=E, run_checks=run_checks)

    if print_output:
        print("Optimal Query Complexity:", np.round(solution['query_complexity'], round_to))
        print("Number of Iterations:", solution['num_iteration'])
        print("Run Time:", np.round(solution['run_time'], round_to), "seconds")
        print()

    return solution

def runSDPForN(getD, getE, n_end, n_start=1, print_output=True, round_to=3, run_checks=3):
    '''
        Parameters:
            getD : function that returns Boolean input D and 
                   takes length of bitstring n as input
            getE : function that returns Boolean output E and
                   takes D as input
            n_end : most bits to run SDP on
            n_start : fewest bits to run SDP on
            print_output : whether to print output of wrapSPDSolver
            round_to : accuracy of query complexity and time (seconds)
            run_checks : whether to run checks on span program
        Returns:
            solutions : dictionary with metadata from call to SDP
                query_complexity : optimal quantum query complexity
                                   of function from D to E for each n
                matrix_solution : array of solution X to SDP for each n
                num_iteration : number of iterations solveSDP ran for for each n
                span_vectors : input vectors of span program for each n
                target_vector : target vector of span program for each n
                run_time : time (seconds) to run whole function for each n
                n_bitstring : n
    '''
    solutions = {}
    for n in range(n_start, n_end + 1):
        D = getD(n=n)
        E = getE(D=D)
        solution = runSDP(D=D, E=E, print_output=print_output, round_to=round_to, run_checks=run_checks)
        for key in solution:
            if key not in solutions:
                solutions[key] = []
            solutions[key] += [solution[key]]
    return solutions
