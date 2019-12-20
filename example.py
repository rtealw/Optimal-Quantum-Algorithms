import quantum_query_optimizer as qqo
import numpy as np

# Example 1
D = ['00', '01', '10', '11']
E = ['0', '1', '1', '1']
qqo.runSDP(D=D, E=E)

# Example 2
solutions = qqo.runSDPForN(getD=qqo.getDAll, getE=qqo.getERandom, n_end=2, n_start=2)

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
