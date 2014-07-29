Work in Progress
===============

 ....



.. code-block:: julia

  using Mopt

  # get a parameter vector
  p = ["a" => 3.1 , "b" => 4.9]
  # define params to use with bounds
  pb= [ "a" => [0,1] , "b" => [0,1] ]

  # get some moments
  # first entry is moment estimate, second is standard deviation
  moms = [
    "alpha" => [ 0.8 , 0.02 ],
    "beta"  => [ 0.8 , 0.02 ],
    "gamma" => [ 0.8 , 0.02 ]
  ]

  # a subset of moments to match
  submoms = ["alpha", "beta"]

  # call objective
  x = Mopt.Testobj(p,moms,submoms)

  # Define an Moment Optimization Problem
  mprob = Mopt.MProb(p,pb,Mopt.Testobj,moms;moments_subset=submoms)

  # show
  mprob

  # step 2: choose an algorithm
  # ----------------------
  algo = Mopt.MAlgoRandom(mprob,opts=["mode"=>"serial","maxiter"=>100])

  # step 3: run estimation
  # ----------------------
  runMopt(algo)

.. py:function:: MProb(ss)
  
  this creates an MProb object
  
  :param str ss: some argument to the function

  
# Efficient Kronecker Product Computation
# ----------------------------------------

.. py:function:: getTensorCoef(ibm,v)

   estimates approximating coefficients

   :param str ibm: Dict{Integer,Array{T,2}} of inverse basis functions
   :param str v: vector of function values

This function is useful for approximation high dimensional functional spaces with basis functions. Suppose we want to approximate an M-dimensional object \eqn{f} that maps R^M into R,  on the tensor product \code{A} of univariate grids \eqn{x_i} of length \eqn{n_i,i=1,\dots,M} each. Denoting the size \eqn{n_i x n_i} basis matrix for dimension \eqn{i} by \eqn{mat_i}, this kronecker product is \code{A <- kronecker(mat_M,kronecker(mat_M-1,...,kronecker(mat_2,mat_1))...)}.\cr
  Approximating coefficients \code{coefs} may be obtained by solving the system \code{A * coefs = y}, where \eqn{y_{i,j,...,k} = f(x_{1,i},x_{2,j},\dots,x_{M,k})} are function values at the grid. Beyond a certain number of dimensions M, or number of data points \eqn{n_i} this is infeasible. Firstly because \code{A} does not fit in memory, and secondly because solving the system is very costly. \code{armakron} implements the efficient deBoor (1979) algorithm to compute \code{coefs}. The kronecker product is never formed, and the calculation involves a series of repeated matrix multiplications. To obtain coefficients, it is most efficient to supply \emph{inverse} basis matrices. This algorithm is orders of magnitude faster than for example setting up \code{A} as a product of sparse matrices (as is the case with spline basis), and to use \code{solve(A,y)} from a sparse library.






