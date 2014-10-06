


# ApproXD.jl: Function Approximation in X Dimensions in Julia 

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


Documentation is started to [being put together on readthedocs](http://approxdjl.readthedocs.org/en/latest/index.html). The best way to get started is currently to look at the unit tests at the moment. sorry about that.


The package contains kronecker products of univariate spline basis functions (tested up to 10D), where the approximating coefficients are computed in an efficient manner according to the algorithm outlined in [1]. 

Next to that the package contains a high performance linear interpolator with a caching mechanism similar to what you get from a [GSL interpolation accelerator](https://www.gnu.org/software/gsl/manual/html_node/Index-Look_002dup-and-Acceleration.html). This supports up to 4D for now.

[1] C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173–182, 1979.

