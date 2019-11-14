# OptQuant
### By R. Teal Witter and Michael T. Czekanski
#### A toolkit to find the optimal quantum query complexity and query optimal quantum algorithm of arbitrary Boolean functions.


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
(Note: We provide example functions in `scripts/handle_input`.)

## Example 1
We consider the Boolean function `OR` on input bitstrings of length 2.
The output is `'1'` if any bit is 1 and `'0'` otherwise.
In this example, we explicitly define both `D` and `E`.
Then we call our function `runSDP` from `scripts/wrap_adm.py`.

```python
D = ['00', '01', '10', '11']
E = ['0', '1', '1', '1']
runSDP(D=D, E=E)
```
The corresponding output should look similar to:
```
n: 2
D: ['00', '01', '10', '11']
E: ['0', '1', '1', '1']
Optimal Query Complexity: 1.414
Number of Iterations: 73
Run Time: 0.067 seconds
```

## Example 2
We again consider `OR` on bitstrings of length 2.
In this example though, we define functions to generate
all bitstrings of length n and evaluate the function `OR` on D.
Then we pass our functions into `runSDPForN` from `scripts/wrap_adm.py`
and specify the starting and ending `n`.
```python
def getDAll(n):
    return [np.binary_repr(i, width=n) for i in range(2**n)]

def getEOR(D):
    return ['1' if '1' in x else '0' for x in D]

runSDPForN(getD=getDAll, getE=getEOR, n_end=2, n_start=2))
```
The corresponding output should look similar to:
```
n: 2
D: ['00', '01', '10', '11']
E: ['0', '1', '1', '1']
Optimal Query Complexity: 1.414
Number of Iterations: 73
Run Time: 0.058 seconds
```

## Semidefinite Program Formulation
We use Ben Reichardt's formulation of the semidefinite program (SDP) for
optimal quantum query complexity (described in `Theorem 6.2`) 
and query optimal span program (`Lemma 6.5`) in
[Span programs and quantum query complexity:
The general adversary bound is nearly tight for every boolean function](https://arxiv.org/pdf/0904.2759.pdf).

## Alternating Direction Method
To solve Reichardt's SDP,
we use Zaiwen Wen, Donald Goldfarb, and Wotao Yin's
`Algorithm 1` described in
[Alternating direction augmented Lagrangian methods for semidefinite programming](http://mpc.zib.de/index.php/MPC/article/viewFile/40/20).