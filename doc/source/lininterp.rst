.. _lininterp:

Linear Interpolator
===================

The ``lininterp`` type is an up to 4-dimensional linear interpolator (and users with higher dimension requirement should note that extending to higher dimensions is relatively easy). 
The main features are

* A caching mechanism: 
	when interpolating ``f(x)`` on a grid of points ``X``, upon being given a point `y` the interpolator needs to find the correct bin or bracket of ``X`` in which ``y`` is contained. This operation can be sped up dramatically if `y` lies in the same bracket as the previous evaluation `z`. The implementation here follows closely `the GSL interpolation accelerator object <https://www.gnu.org/software/gsl/manual/html_node/Index-Look_002dup-and-Acceleration.html>`_
* The possbility to evaluate mutiple functions on the same point in a single pass:
	if you have multiple function ``y1,...,yk``, you can contruct the ``lininterp`` object by supplying one array of function values for each of them, and select the appropriate index at the time you want to obtain function ``j``, say.

Construction
------------



Methods
-------


.. function::getValue(l,point)

standard evaluator: evaluates all functions stored on ``l`` at ``point``.
``point`` is a ``Vector{Float64}``.

.. function::getValue!(y,l,point,which)






 py:function:: lininterp(v,g)

 Multidimensional Linear Interpolator with cache accelerator. Can take multiple functions to be approximated 
 on the same point. 
 ``v`` is an array of function values (i.e. ``Array{Float64}``), or an array with several function value arrays (``Array{Array{Float64}}``).
 ``gs`` is an Array of one-dimensional grids. 
 
 The dimensions of ``gs`` and ``v`` must correspond, in the sense that ``size(v,i) == length(gs[i])``, i.e. dimension `i` of the function array must be as long as the grid for this dimension.
