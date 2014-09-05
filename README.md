


# ApproXD: Function Approximation in X Dimensions in Julia 

[![Build Status](https://travis-ci.org/floswald/ApproXD.jl.png?branch=master)](https://travis-ci.org/floswald/ApproXD.jl)

<!-- [![Coverage Status](https://coveralls.io/repos/floswald/BSplines.jl/badge.png)](https://coveralls.io/r/floswald/BSplines.jl) -->


[![2D approximation](https://dl.dropboxusercontent.com/u/109115/BSplines.jl/approx.png)]()

This package implements high-dimensional approximation of real-valued functions in Julia:

```julia
x = [1:10.0]
y = [1:10.0]
z = [1:10.0]
fun(x,y,z) = x^2 + y^(0.5) + (x-z)^2

# what is fun(4.5,pi,8.01) ?
```

The package aims in particular at a fast solution to this problem. The focus is to get an interpolated point from a high-dimensional object.

At the moment it contains kronecker products of univariate spline basis functions (tested up to 10D), where the approximating coefficients are computed in an efficient manner according to the algorithm outlined in [1]. 

Next to that the package contains a high performance linear interpolator with a cashing mechanism similar to what you get from a [GSL interpolation accelerator](https://www.gnu.org/software/gsl/manual/html_node/Index-Look_002dup-and-Acceleration.html). This supports only 3D up to now, but extension is in progress.

[1] C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173–182, 1979.

