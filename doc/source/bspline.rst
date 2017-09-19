
.. _bspline-label:

BSpline
=======

Provides a simple Bspline basis function type with two constructors. One with equally spaced knots, and one with user set knots. The knot vectors are *extended* knot vectors, i.e. for a degree 3 spline with interior knots (3,4,5), we would have knot vector (3,3,3,3,4,5,5,5,5).

Construction
------------

.. function:: BSpline(nKnots,deg,lb,ub)

	``nKnots``: number of knots

	``deg``: Spline degree

	``lb``: lower bound of spline space

	``ub``: upper bound

	This constructor allocates ``2*deg + nKnots`` equally spaced knots. 


.. function:: BSpline(knots,deg)

	``knots``: knot vector

	``deg``: Spline degree

	Takes user supplied knot vector. Do **not** supply the extended vector.

