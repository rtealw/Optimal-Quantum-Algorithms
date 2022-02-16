import numpy as np
import quantum_query_optimizer as qqo

def test_getDAll3():
    n=3
    D = qqo.getDAll(n)

    truth = ['000', '001', '010', '011', '100', '101', '110', '111']
    assert(sum([x in truth for x in D]) == len(truth))

def test_getDWorstOR():
    n=3
    D = qqo.getDWorstOR(n)

    truth = ['000', '001', '010', '100']
    assert(sum([x in truth for x in D]) == len(truth) and len(D) == len(set(D)))


def test_getEOR():
    n=3
    D = qqo.getDWorstOR(n)
    E = qqo.getEOR(D)

    truth = ['0', '1', '1', '1']
    assert(sum([E[i] == truth[i] for i in range(n+1)]) == len(E) and len(E) == len(truth)) 
