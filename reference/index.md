# Package index

## Discrete calculus

Divided differencing, discrete differentiation, and discrete
integration.

- [`divided_diff()`](https://glmgen.github.io/dspline/reference/divided_diff.md)
  : Divided differencing
- [`discrete_deriv()`](https://glmgen.github.io/dspline/reference/discrete_deriv.md)
  : Discrete differentiation
- [`discrete_integ()`](https://glmgen.github.io/dspline/reference/discrete_integ.md)
  : Discrete integration

## Matrix multiplication

Multiplication by discrete derivative and falling factorial basis
matrices.

- [`d_mat_mult()`](https://glmgen.github.io/dspline/reference/d_mat_mult.md)
  : Multiply by D matrix
- [`b_mat_mult()`](https://glmgen.github.io/dspline/reference/b_mat_mult.md)
  : Multiply by B matrix
- [`h_mat_mult()`](https://glmgen.github.io/dspline/reference/h_mat_mult.md)
  : Multiply by H matrix

## Matrix construction

Construction of discrete derivative and discrete spline basis matrices.

- [`d_mat()`](https://glmgen.github.io/dspline/reference/d_mat.md) :
  Construct D matrix
- [`b_mat()`](https://glmgen.github.io/dspline/reference/b_mat.md) :
  Construct B matrix
- [`h_mat()`](https://glmgen.github.io/dspline/reference/h_mat.md) :
  Construct H matrix
- [`n_mat()`](https://glmgen.github.io/dspline/reference/n_mat.md) :
  Construct N matrix

## Basis evaluation

Evaluation of falling factorial and discrete B-spline basis functions.

- [`h_eval()`](https://glmgen.github.io/dspline/reference/h_eval.md) :
  Evaluate H basis
- [`n_eval()`](https://glmgen.github.io/dspline/reference/n_eval.md) :
  Evaluate N basis

## Interpolation

Interpolation within the “canonical” space of discrete splines.

- [`dspline_interp()`](https://glmgen.github.io/dspline/reference/dspline_interp.md)
  : Discrete spline interpolation

## Projection

Least squares projection onto “custom” spaces of discrete splines.

- [`dspline_solve()`](https://glmgen.github.io/dspline/reference/dspline_solve.md)
  : Discrete spline projection

## In-place computations

Divided differencing and matrix multiplication using in-place
operations.

- [`.divided_diff()`](https://glmgen.github.io/dspline/reference/dot_functions.md)
  [`.b_mat_mult()`](https://glmgen.github.io/dspline/reference/dot_functions.md)
  [`.h_mat_mult()`](https://glmgen.github.io/dspline/reference/dot_functions.md)
  : In-place computations
