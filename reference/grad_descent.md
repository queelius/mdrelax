# Gradient descent optimizer with numerical gradients

Minimizes a function `f` using gradient descent with numerical gradient
approximation. Supports constraints via a support function.

## Usage

``` r
grad_descent(
  f,
  x0,
  sup = function(x) TRUE,
  lr = 0.01,
  eps = 1e-08,
  max_iter = 10000L,
  debug = FALSE,
  h = 1e-06
)
```

## Arguments

- f:

  Function to minimize. Takes a vector and returns a scalar.

- x0:

  Initial parameter values.

- sup:

  Support function. Returns TRUE if parameters are valid.

- lr:

  Learning rate (step size).

- eps:

  Convergence tolerance for gradient norm.

- max_iter:

  Maximum number of iterations.

- debug:

  If TRUE, print debugging information.

- h:

  Step size for numerical gradient approximation.

## Value

A list with components:

- param:

  Final parameter values

- converged:

  TRUE if converged within max_iter

- iter:

  Number of iterations performed

- value:

  Final function value
