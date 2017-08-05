# Adaptive Quadrature with Bag of Task

Adaptive quadrature is a recursive algoritm that approximate the intehral of a given function f(x). We implement a parallelized MPI version of the algorithm using the Bag of Tasks pattern. In the program, the function to be computed is

```
cosh(arg)*cosh(arg)*cosh(arg)*cosh(arg)
```

from a=0 and b=5, with epsilon 1e-3. We use one farmer and three workers. The final result is 7583461.801486.

## Reference
[PPLS course page](http://www.inf.ed.ac.uk/teaching/courses/ppls/)
