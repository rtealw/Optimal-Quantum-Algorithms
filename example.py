import sys
sys.path.append('scripts/')
import quantum_query_optimizer as qqo
import boolean_functions as bf

## Example 1
#D = ['00', '01', '10', '11']
#E = ['0', '1', '1', '1']
#qqo.runSDP(D=D, E=E)

## Example 2
#qqo.runSDPForN(getD=bf.getDAll, getE=bf.getEOR, n_end=2, n_start=2)

# solutions = qqo.runSDP(D=['00', '01', '10', '11'], E=['0', '0', '0', '1'])
solutions = qqo.runSDPForN(getD=bf.getDAll, getE=bf.getERandom, n_end=4, n_start=2)