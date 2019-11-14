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
Second, create one function that generates the set D for arbitrary bitstring length n
and another function that generates the set E from D according to f.

## Example 1

```python
D = ['00', '01', '10', '11']
E = ['0', '1', '1', '1']
runSDP(D=D, E=E)
```

```
n: 2
D: ['00', '01', '10', '11']
E: ['0', '1', '1', '1']
Optimal Query Complexity: 1.414
Number of Iterations: 73
Run Time: 0.067 seconds
```

## Example 2

```python
def getDAll(n):
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def getEOR(D):
    return ['1' if '1' in x else '0' for x in D]

runSDPIterations(iterations=2, getD=getDAll, getE=getEOR, start=2)
```

```
n: 2
D: ['00', '01', '10', '11']
E: ['0', '1', '1', '1']
Optimal Query Complexity: 1.414
Number of Iterations: 73
Run Time: 0.058 seconds
```


## Semidefinite Program Formulation

## Alternating Direction Method