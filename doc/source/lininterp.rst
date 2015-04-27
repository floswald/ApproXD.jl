.. _lininterp:

Linear Interpolator
===================

The ``lininterp`` type is an up to 4-dimensional linear interpolator (and users with higher dimension requirement should note that extending to higher dimensions is relatively easy). It's a linear interpolator with cache accelerator. Can take multiple functions to be approximated on the same point.
The main features are

* A caching mechanism: 
	when interpolating ``f(x)`` on a grid of points ``X``, upon being given a point `y` the interpolator needs to find the correct bin or bracket of ``X`` in which ``y`` is contained. This operation can be sped up dramatically if `y` lies in the same bracket as the previous evaluation `z`. The implementation here follows closely `the GSL interpolation accelerator object <https://www.gnu.org/software/gsl/manual/html_node/Index-Look_002dup-and-Acceleration.html>`_
* The possbility to evaluate mutiple functions on the same point in a single pass:
	if you have multiple functions ``f1,...,fk``, you can contruct the ``lininterp`` object by supplying one array of function values for each of them. Depending on the method you use, ``lininterp`` will return the approximation for all functions at ``x``, or only a specific one. 

What ``lininterp`` **cannot** do:

* Multidimensional extrapolation: I am interpolating the supplied function on the convex hull of grid points. in 2D, think of a rectangle that is formed by 4 points ``(x1,y1), (x2,y1),(x1,y2),(x2,y2)``, and some point ``(x,y)`` that lies within. I approximate this point as a convex combination of the four *vertices*, where the weights are the position of ``x`` in ``[x1,x2]`` and that of ``y`` in ``[y1,y2]``. This means I cannot know the value for ``x>x2``, for example. 

Constructors
------------

.. function:: lininterp(v,g)

	:param v: 	an ``Array{Float64,nD}`` of ``n`` function evaluations. ``n=n1*n2*...nD`` for ``nD`` dimensions. 

	:param g: 	an ``Array{Vector{Float64}}``, i.e. a collection of one-dimensional grids. can be irregularly spaced, but must be sorted ascendingly.

	For example, 

	.. code-block:: julia

		v = rand(5,6)				# function values
		g = Array{Float64,1}[]	 	# grid container
		push!(g,linspace(0,1,5),linspace(-1.1,-0.1,6))  # grids
		L = lininterp(v,g)   		# interpolator 

.. function:: lininterp(v1,v2,g)

	:param v1: an ``Array{Float64,nD}`` of ``n`` function evaluations for function ``f1``. ``n=n1*n2*...nD`` for ``nD`` dimensions.

	:param v2: an ``Array{Float64,nD}`` of ``n`` function evaluations for function ``f2``. ``n=n1*n2*...nD`` for ``nD`` dimensions.

	:param g: 	an ``Array{Vector{Float64}}``, i.e. a collection of one-dimensional grids. can be irregularly spaced, but must be sorted ascendingly.

	Note that ``v1,v2`` are defined on the same grid ``g``.

	For example, 

	.. code-block:: julia

		v1 = rand(5,6)				# function values 1
		v2 = rand(5,6)				# function values 2
		g = Array{Float64,1}[]	 	# grid container
		push!(g,linspace(0,1,5),linspace(-1.1,-0.1,6))  # grids
		L = lininterp(v1,v2,g)   		# interpolator 

.. function:: lininterp(v,g)

	:param v: a collection of multiple function evaluations defined on the same grid, i.e. an ``Array{Array{Float64,nD}}``

	:param g: 	an ``Array{Vector{Float64}}``, i.e. a collection of one-dimensional grids. can be irregularly spaced, but must be sorted ascendingly.

	Note that all ``v`` are defined on the same grid ``g``.

	For example, 

	.. code-block:: julia

		v1 = rand(3,5,6)			# function values 1
		v2 = rand(3,5,6)			# function values 2
		v3 = rand(3,5,6)			# function values 3
		v = Array{Float64}[]		# values container
		g = Array{Float64,1}[]	 	# grid container
		push!(v,v1,v2,v3) 			# grids
		push!(g,linspace(0,1,3),linspace(-1.1,-0.1,5),linspace(2,4.1,6))  # grids
		L = lininterp(v,g)   		# interpolator 


Methods
-------

.. function:: getValue(L,point)

	:param L: 	a ``lininterp`` object

	:param point: 	coordinates of a point, a vector

	returns the approximation for all functions stored in ``L`` at ``point``.

.. function:: getValue!(y,L,point,which)

	:param y: 	prealloacted return value. vector which ``length(which)``

	:param L: 	a ``lininterp`` object

	:param point: coordinates of a point (a vector)

	:param which: integer index of which function to be evaluated

	if there are multiple functions stored in ``l``, defined on the same grids, you can select ``which`` one will be evaluated on ``point`` by setting ``which``.







