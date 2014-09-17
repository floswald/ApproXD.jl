.. ApproXD documentation master file, created by
   sphinx-quickstart on Wed Sep 17 11:03:28 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ApproXD's documentation!
===================================

ApproXD.jl is a Julia package for function approximation in high dimensional space. The problem at hand is that a real-valued function

.. math::
	f(x), x \in R^d

is known on a finite set of points :math:`X \subset R^d`, but not in between. The task is to provide an interpolation (*connecting the dots*) that makes evaluation off the grid :math:`X` possible. There are **two supported ways to make the interpolation**, Tensor Product of Spline Basis functions and Linear Interpolator.

Contents:

.. toctree::
   :maxdepth: 2

   lininterp.rst
   tensor.rst
   bspline.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

