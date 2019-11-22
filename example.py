import sys
sys.path.append('scripts/')

import quantum_query_optimizer as qqo
import matplotlib.pyplot as plt

## Example 1
#D = ['00', '01', '10', '11']
#E = ['0', '1', '1', '1']
#qqo.runSDP(D=D, E=E)

## Example 2
#qqo.runSDPForN(getD=bf.getDAll, getE=bf.getEOR, n_end=2, n_start=2)

# solutions = qqo.runSDP(D=['00', '01', '10', '11'], E=['0', '0', '0', '1'])

# Generate paper figures
all_solutions = qqo.runSDPForN(getD=qqo.getDAll, getE=qqo.getEOR, n_end=6, n_start=1)
qqo.visualizeRuntime(all_solutions, title="Runtime of OR by Input Size for All Inputs", filename="figures/or_all_runtime.png")
qqo.visualizeComplexity(all_solutions, title="Complexity of OR by Input Size for All Inputs", filename="figures/or_all_complexity.png")

worst_solutions = qqo.runSDPForN(getD=qqo.getDWorstOR, getE=qqo.getEOR, n_end=20, n_start=1)
qqo.visualizeRuntime(worst_solutions, title="Runtime of OR by Input Size for Worse Inputs", filename="figures/or_worst_runtime.png")
qqo.visualizeComplexity(worst_solutions, title="Complexity of OR by Input Size for Worse Inputs", filename="figures/or_worst_complexity.png")
