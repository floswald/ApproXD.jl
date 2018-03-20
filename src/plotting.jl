

using PyPlot



	# 2D approximation example
function plot2D()

	ndims = 2

	doplot = false

	# bounds
	lb = [-2.1,-2.5]
	ub = [2.2,2.6]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [11,12]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => collect(linspace(lb[i],ub[i],npoints[i])) for i=1:ndims)

	# set up BSplines
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

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
	function f(x,y) 
		sin(sqrt(x^2+y^2))
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j) for i in points[1], j in points[2]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 1e-2


	 # make a plot
	 # ===========

	new_npoints = [17,16]
	# new_npoints = npoints
	new_points = Dict(i => linspace(lb[i],ub[i],new_npoints[i]) for i=1:ndims)
	# set of new basis functions
	nd = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		nd[i] = full(getBasis(new_points[i],bsp[i]))
	end

	pred = BSplines.evalTensor2(nd[2],nd[1],mycoef)
	
	BSplines.figure()
	BSplines.subplot(121,projection="3d")
	BSplines.mesh(points[2],points[1],y)	
	BSplines.title("true values on grid")

	BSplines.subplot(122,projection="3d")
	pp = reshape(pred,new_npoints[1],new_npoints[2])
	BSplines.mesh(new_points[2],new_points[1],pp)	
	BSplines.title("Approximation off grid")

	BSplines.suptitle("2D TEST")
end


function plot3D()	
	# 3D approximation example

	ndims = 3
	doplot = false

	# bounds
	lb = [-2.1,-2.5,0.9]
	ub = [2.2,2.6,5.0]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [12,16,13]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,3,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => linspace(lb[i],ub[i],npoints[i]) for i=1:ndims)

	# set up BSplines
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

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
	function f(x,y,z) 
		sin(sqrt(x^2+y^2)) - z^2
	end

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j,k) for i in points[1], j in points[2], k in points[3]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)


		 # make a plot
		 # ===========

		new_npoints = [17,16,5]
		# new_npoints = npoints
		new_points = Dict(i => linspace(lb[i],ub[i],new_npoints[i]) for i=1:ndims)
		# set of new basis functions
		nd = Dict{Integer,Array{Float64,2}}()
		for i=1:ndims
			nd[i] = full(getBasis(new_points[i],bsp[i]))
		end

		pred = BSplines.evalTensor3(nd[3],nd[2],nd[1],mycoef)
		
		BSplines.figure()
		BSplines.subplot(221,projection="3d")
		BSplines.mesh(points[2],points[1],y[:,:,1])	
		BSplines.title("true values at lowest 3D state")

		BSplines.subplot(222,projection="3d")
		pp = reshape(pred,new_npoints[1],new_npoints[2],new_npoints[3])
		BSplines.mesh(new_points[2],new_points[1],pp[:,:,1])	
		BSplines.title("Approximation")

		BSplines.subplot(223,projection="3d")
		BSplines.mesh(points[2],points[1],y[:,:,npoints[3]])	
		BSplines.title("true values at highest 3D state")

		BSplines.subplot(224,projection="3d")
		pp = reshape(pred,new_npoints[1],new_npoints[2],new_npoints[3])
		BSplines.mesh(new_points[2],new_points[1],pp[:,:,new_npoints[3]])	
		BSplines.title("Approximation")
		BSplines.suptitle("3D TEST")

end


function plot4D()	
	ndims = 4

	# bounds
	lb = [-1.1,-1.5,-0.9,-1.0]
	ub = [1.2,1.6,0.9,1]

	# number of eval points and basis functions:
	# we require square basis matrices!
	npoints = [13,13,13,13]

	# number of basis funcs
	nbasis = npoints

	# splien degrees
	degs = [3,3,3,3]

	# implies a number of knots for each spline
	# remember the restriction that nknots == ncoefs
	nknots = Dict(i => nbasis[i] - degs[i] + 1 for i=1:ndims)

	# eval points
	points = Dict(i => linspace(lb[i],ub[i],npoints[i]) for i=1:ndims)

	# set up BSplines
	bsp = Dict(i => BSpline(nknots[i],degs[i],lb[i],ub[i]) for i=1:ndims)

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

	# compute function values so that
	# k is fastest varying index
	y = Float64[f(i,j,k,w) for i in points[1], j in points[2], k in points[3], w in points[4]]

	yvec = y[:]

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	 # testing tolerance
	 tol = 5e-3


		 # make a plot
		 # ===========

		new_npoints = [17,16,5,7]
		# new_npoints = npoints
		new_points = Dict(i => linspace(lb[i],ub[i],new_npoints[i]) for i=1:ndims)
		# set of new basis functions
		nd = Dict{Integer,Array{Float64,2}}()
		for i=1:ndims
			nd[i] = full(getBasis(new_points[i],bsp[i]))
		end

		# predict everywhere
		pred = BSplines.evalTensor4(nd[4],nd[3],nd[2],nd[1],mycoef)
		
		BSplines.figure()
		BSplines.subplot(221,projection="3d")
		BSplines.mesh(points[2],points[1],y[:,:,1,1])	
		BSplines.title("true values at lowest state (:,:,1,1)")

		BSplines.subplot(222,projection="3d")
		pp = reshape(pred,new_npoints[1],new_npoints[2],new_npoints[3],new_npoints[4])
		BSplines.mesh(new_points[2],new_points[1],pp[:,:,1,1])	
		BSplines.title("Approximation")

		BSplines.subplot(223,projection="3d")
		BSplines.mesh(points[2],points[1],y[:,:,npoints[3],npoints[4]])	
		BSplines.title("true values at state (:,:,end,end)")

		BSplines.subplot(224,projection="3d")
		pp = reshape(pred,new_npoints[1],new_npoints[2],new_npoints[3],new_npoints[4])
		BSplines.mesh(new_points[2],new_points[1],pp[:,:,new_npoints[3],new_npoints[4]])	
		BSplines.title("Approximation")
		BSplines.suptitle("4D TEST")
end