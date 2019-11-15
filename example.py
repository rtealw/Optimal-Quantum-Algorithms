import sys
sys.path.insert(1, 'scripts/')

from boolean_functions import *
import quantum_query_optimizer as qqo

## Example 1
#D = ['00', '01', '10', '11']
#E = ['0', '1', '1', '1']
#qqo.runSDP(D=D, E=E)

## Example 2
#qqo.runSDPForN(getD=getDAll, getE=getEOR, n_end=2, n_start=2)

solutions = qqo.runSDPForN(getD=getDAll, getE=getEOR, n_end=5, n_start=1)