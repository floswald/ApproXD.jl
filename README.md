# ApproXD

[![Build Status](https://travis-ci.com/floswald/ApproXD.jl.svg?branch=master)](https://travis-ci.com/floswald/ApproXD.jl)
[![Codecov](https://codecov.io/gh/floswald/ApproXD.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/floswald/ApproXD.jl)


* This package implements bspline and linear interpolation in julia
* For most purposes, [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) will be preferrable to this package.
* However, there are some features which are available here, and not there.
    - the method `getTensorCoef` is a very efficient algorithm to compute approximating coefficients from a tensor product of basis matrices. it is efficient because it never forms the tensor product.
    - the package allows low-level access to objects such as spline knot vectors. Suppose you want to have a knot vector with a [knot multiplicity](https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/B-spline/bspline-mod-knot.html) in the interior knot span to approximate a kink. For example,
    ```julia
    knots = vcat(lb,-0.5,0,0,0.5,ub)
    ```
    is a valid knot vector.
* Documentation is non-existent. Please look at the tests. Sorry.

