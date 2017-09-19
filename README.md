


# ApproXD.jl: Function Approximation in X Dimensions in Julia 

[![Build Status](https://travis-ci.org/floswald/ApproXD.jl.png?branch=master)](https://travis-ci.org/floswald/ApproXD.jl)

[![Build status](https://ci.appveyor.com/api/projects/status/p4tr6m340xa1r9a6?svg=true)](https://ci.appveyor.com/project/floswald/approxd-jl)

[![Coverage Status](https://coveralls.io/repos/floswald/ApproXD.jl/badge.png)](https://coveralls.io/r/floswald/ApproXD.jl)



This package implements high-dimensional approximation of real-valued functions in Julia:

```julia
x = [1:10.0]
y = [1:10.0]
z = [1:10.0]
fun(x,y,z) = x^2 + y^(0.5) + (x-z)^2

# what is fun(4.5,pi,8.01) ?
```

