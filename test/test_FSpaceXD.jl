
module test_fspace
using ApproXD
using Test



@testset "testing FSpaceXD" begin
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
	points = Dict(i => collect(range(lb[i],stop = ub[i], length = npoints[i])) for i=1:ndims)

	# set up ApproXD
	bsp = Dict{Integer,BSpline}()
	for i in 1:ndims
		bsp[i] = BSpline(nknots[i],degs[i],lb[i],ub[i])
	end

	# set of basis functions
	d = Dict{Integer,Array{Float64,2}}()
	for i=1:ndims
		d[i] = Array(getBasis(points[i],bsp[i]))
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

	# get coefs using the function
	mycoef = getTensorCoef(id,yvec)

	# setup the FSpace
	fx = FSpaceXD(ndims,mycoef,bsp)
	
	rval1 = lb[1] + 0.3 
	rval2 = lb[2] + 0.23
	rval3 = lb[3] + 0.111
	rval4 = lb[4] + 0.099
	# println("approx value = $(getValue(fx,[rval1,rval2,rval3,rval4]))")
	# println("true value = $(f(rval1,rval2,rval3,rval4))")
	@test isapprox(getValue(fx,[rval1,rval2,rval3,rval4]), f(rval1,rval2,rval3,rval4),atol=3e-3)


end


end