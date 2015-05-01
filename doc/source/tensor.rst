Tensor Product of Spline Basis functions
========================================


The main user interface to approximate a function is called ``FSpaceXD``. There are several lower level functions which may be of interest to users, and they will be documented below.

Construction of ``FSpaceXD``
----------------------------

.. function:: FSpaceXD(ndim,coeff,bs)

	``ndim``: number of dimensions

	``coeff``: vector of approximating coefficients

	``bs``: ``Dict`` of length ``ndim`` with a spline basis function for each dimension. 


Methods for ``FSpaceXD``
------------------------

.. function:: getValue(fspace,x)
    :noindex:

    ``fspace``: an ``FSpaceXD`` object

    ``x``: a vector of points, ``length(x)=ndim``



``FSpaceXD`` Useage Example
----------------------------

the example uses the :ref:`bspline-label` type.

.. code-block:: julia

	# set ndims
	ndims = 4

	# bounds
	lb = [-1.1,-1.5,-0.9,-1.0]
	ub = [1.2,1.6,0.9,1]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [10,11,9,8]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,1,3,2]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = [i => nbasis[i] - degs[i] + 1 for i=1:ndims]

	# eval points
	points = [i => linspace(lb[i],ub[i],npoints[i]) for i=1:ndims]

	# set up ApproXD
	bsp = Dict{Integer,BSpline}()
	for i in 1:ndims
		bsp[i] = BSpline(nknots[i],degs[i],lb[i],ub[i])
	end

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = full(getBasis(points[i],bsp[i]))
	end

	# set of INVERSE basis functions
	id = Dict{Integer,Array{Float64,2}}()
	for k in collect(keys(d))
		id[k] = inv(d[k])
	end

	#  get a function
	function f(x,y,z,w) 
		sin(sqrt(x^2+y^2)) + (z-w)^3
	end

	y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

	yvec = y[:]

	# get coefs using the Tensor approximator function
	mycoef = getTensorCoef(id,yvec)

	# setup the FSpace
	fx = FSpaceXD(ndims,mycoef,bsp)
	
	rval1 = lb[1] + 0.3 
	rval2 = lb[2] + 0.23
	rval3 = lb[3] + 0.111
	rval4 = lb[4] + 0.099
	println("approx value = $(getValue(fx,[rval1,rval2,rval3,rval4]))")
	println("true value = $(f(rval1,rval2,rval3,rval4))")
	




Outline of the de Boor Algorithm
--------------------------------

The main features of this implementation are 

1. Avoids allocation of large tensor product
2. Exploits sparsity of spline basis functions

For a point of interest :math:`x \in R^d`, form the tensor product of univariate Bsplines of arbitrary degree ``k`` on each grid dimension and estimate a coefficient vector :math:`c` by solving the following system:

.. math::
	c &= \Phi(X_1,X_2,\dots,X_d)^{-1} f(X) \\
	X_1 &= \{x_{11},x_{12},\dots,x_{1n(1)} \}\\
	\dots \\
	X_d &= \{x_{d1},x_{d2},\dots,x_{dn(d)} \} \\
	\Phi(X_1,X_2,\dots,X_d) &= B_1(X_1) \otimes B_2(X_2) \otimes \dots \otimes B_d(X_d) 

Here the :math:`B`'s are spline basis functions of arbitrary degree, and they are evaluated on the entire grid in each dimension. **An important restriction** in this setup is that the resulting Bspline basis function matrices must be square, i.e. there must be as many evaluation points as there are resulting basis functions (given ``n`` evaluation points and a selected spline degree ``k``, this results in a predetermined number of spline knots and therefore basis functions.)

Once the coefficients are obtained, the function value is obtained via

.. math::
	\hat{f}(x_1,x_2,\dots,x_d) = c \times \Phi(x_1,x_2,\dots,x_d)

The problem here is that :math:`\Phi` quickly becomes very large; in fact so large that it becomes infeasible to even allocate it in memory, let alone solve the system. The approach implmented here never forms :math:`\Phi` and uses a very efficient approach to to solve the system. The algorithm is devised and described in 

	C. De Boor. Efficient computer manipulation of tensor products. ACM Transactions on Mathematical Software (TOMS), 5(2):173–182, 1979.





