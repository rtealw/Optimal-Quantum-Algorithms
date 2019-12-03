import numpy as np
from quantum_query_optimizer import wrapSDPSolver
from termcolor import cprint 

def testSDPSolver(n_end=5, accuracy = 2):
    '''
        Parameters:
            n_end : highest number of bit strings
            accuracy : place to round expected and received to
        Returns:
           (printed output) : test summary
    '''
    all_passed = True
    for n in range(1, n_end + 1):
        print("Testing SDP solver on OR for n = {}...".format(n), end=" ")
        D = [np.binary_repr(i, width=n) for i in range(2**n)]
        E = ['1' if '1' in x else '0' for x in D]
        solution = wrapSDPSolver(D=D, E=E, run_checks=True)
        expected = np.sqrt(n)

        # check optimization results
        if round(solution['query_complexity'], accuracy) == round(expected, accuracy):
            cprint("SDP solution passed :)", "green")
        else:
            all_passed = False
            cprint("SDP solution failed :(", "red")
            cprint("Expected: {}".format(expected), "green")
            cprint("Received: {}".format(received), "red")

    if all_passed:
        cprint("\n All tests passed :)", "green")
    else:
        cprint("\n Tests failed :(", "red")

if __name__ == '__main__':
    testSDPSolver()
