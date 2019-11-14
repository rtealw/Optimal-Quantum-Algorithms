# OptQuant
### By R. Teal Witter and Michael T. Czekanski

#### A toolkit to find the optimal quantum query complexity and query optimal quantum algorithm for given Boolean functions.

Consider a function f that maps from D to E where D is a subset of bitstrings
of length n and E is the set of single bit outputs.
In the query model, an algorithm looks at the bits of the input string x in D
as few times as possible before correctly determing f(x).
Given f, our program finds the optimal query complexity of a quantum algorithm
that evaluates f and a span program (i.e. quantum algorithm) that meets
this query complexity.

There are two ways to run our program.
First, explicitly specify the sets D and E.
Second, create functions that generate the set D and evaluate f
for arbitrary bitstring size n.

## Example 1
```python
import handle_input from scripts
D = ['000', '001', '010', '100']
E = ['0', '1', '1', '1']
runSDP(D=D, E=E)
```

## Example 2

## Semidefinite Program Formulation

## Alternating Direction Method