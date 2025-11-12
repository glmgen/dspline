# Divided differencing

Computes the divided difference of a function (or vector of function
evaluations) with respect to given centers.

## Usage

``` r
divided_diff(f, z)
```

## Arguments

- f:

  Function, or vector of function evaluations at the centers.

- z:

  Centers for the divided difference calculation.

## Value

Divided difference of `f` with respect to centers `z`.

## Details

The divided difference of a function \\f\\ with respect to centers
\\z_1, \ldots, z\_{k+1}\\ is defined recursively as: \$\$
f\[z_1,\ldots,z\_{k+1}\] = \displaystyle
\frac{f\[z_2,\ldots,z\_{k+1}\] - f\[z_1,\ldots,z_k\]}{z\_{k+1}-z_1},
\$\$ with base case \\f\[z_1\] = f(z_1)\\ (that is, divided differencing
with respect to a single point reduces to function evaluation).

A notable special case is when the centers are evenly-spaced, say, \\z_i
= z+ih\\, \\i=0,\ldots,k\\ for some spacing \\h\>0\\, in which case the
divided difference becomes a (scaled) forward difference, or
equivalently a (scaled) backward difference, \$\$ k! \cdot
f\[z,\ldots,z+kh\] = \displaystyle \frac{1}{h^k} (F^k_h f)(z) =
\frac{1}{h^k} (B^k_h f)(z+kh), \$\$ where we use \\F^k_h\\ and \\B^k_v\\
to denote the forward and backward difference operators, respectively,
of order \\k\\ and with spacing \\h\\.

## Examples

``` r
f = runif(4)
z = runif(4)
divided_diff(f[1], z[1])
#> [1] 0.4790245
f[1]
#> [1] 0.4790245
divided_diff(f[1:2], z[1:2])
#> [1] -1.281506
(f[1]-f[2])/(z[1]-z[2])
#> [1] -1.281506
divided_diff(f[1:3], z[1:3])
#> [1] 3.748378
((f[1]-f[2])/(z[1]-z[2]) - (f[2]-f[3])/(z[2]-z[3])) / (z[1]-z[3]) 
#> [1] 3.748378
divided_diff(f, 1:4)
#> [1] -0.05887262
diff(f, diff = 3) / factorial(3)
#> [1] -0.05887262
```
