

mutable struct Lininterp

	d          :: Array{Int}    # number of points in each dimension
	n          :: Int 	        # number  of dims
	nfunc      :: Int 	        # number  of different functions to evaluate on the same point
	ifunc      :: Vector{Int} 	# indices of which funcs to evaluate
	hascache   :: Bool 	        # has active cache
	cache      :: Array{Int}	# current index of lower bound of bracket containing x
	infs       :: Array{Float64}	    # current value of inf: lower bound of bracket in xgrid
	sups       :: Array{Float64}	   	# current value of sup: upper bound of bracket in xgrid
	z          :: Array{Float64}    	# current position in bracket: z = (x-x[inf]) / (x[sup] - x[inf])
	hits       :: Int 			# count of cache hits
	miss       :: Int 			# count of cache misses
	grids      :: Array{Vector{Float64}} 	# grids for each dimension
	vals       :: Array{Array{Float64}} 	# arrays with function values on the grid.
	vertex     :: Array{Array{Float64}}  # arrays with current vertices

	# constructor for multiple functions to approximate
	function Lininterp(v::Array{Array{Float64}},g::Array{Vector{Float64}})
		d = [size(v[1])...]
		if !all(map(z->size(z)==size(v[1]),v))
			throw(ArgumentError("all arrays in v must have the same shape"))
		end
		n = length(d)
		nfunc = length(v)
		ifunc = collect(1:nfunc)
		if length(d) != length(g)
			throw(ArgumentError("v must have as many dims as g has grids"))
		end
		for i in 1:length(g)
			if d[i] != length(g[i])
				throw(ArgumentError("size(v,$i) must have the same length as g[$i]"))
			end
			if !issorted(g[i])
				throw(ArgumentError("g[$i] must be sorted!"))
			end
		end
		if n > 4
			throw(ArgumentError("currently only up to 4D implemented"))
		end
		cache = ones(Int,n)
		vertex = Array{Float64}[]
		for i=1:nfunc
			push!(vertex,zeros(2^n))
		end
		return new(d,n,nfunc,ifunc,false,cache,zeros(n),zeros(n),zeros(n),0,0,g,v,vertex)
	end

end

# additional outer constructors
# constructor for single function
function Lininterp(v::Array{Float64},g::Array{Vector{Float64}})
	v1 = Array{Float64}[]
	push!(v1,v)
	Lininterp(v1,g)
end
# constructor for two functions
function Lininterp(v1::Array{Float64},v2::Array{Float64},g::Array{Vector{Float64}})
	v = Array{Float64}[]
	push!(v,v1)
	push!(v,v2)
	Lininterp(v,g)
end


# general utility methods
# =======================

getDims(l::Lininterp) = l.d 
getGrids(l::Lininterp) =  l.grids 
function resetCache!(l::Lininterp) 
	fill!(l.cache,0)
	l.hits = 0
	l.miss = 0
	l.hascache=false
	return nothing
end
function getCache(l::Lininterp)
	l.cache
end
function getCache(l::Lininterp,i::Int)
	l.cache[i]
end
function getCachedVal(l::Lininterp,i::Int)
	l.grids[i][l.cache[i]]
end
function getNextCachedVal(l::Lininterp,i::Int)
	l.grids[i][l.cache[i] + 1]
end

function setValue!(l::Lininterp,which::Int,v::Array)
	if size(v) != l.d
		throw(ArgumentError("size(v) must be equal to size d"))
	end
	vals[which] = v
	return nothing
end

function setGrid!(l::Lininterp,i::Int,g::Vector{Float64})
	if length(g) != l.d[i]
		throw(ArgumentError("new grid g incompatible with dim $i = $(l.d[i])"))
	end
	if !issorted(g)
		throw(ArgumentError("g must be sorted!"))
	end
	l.grids[i] = g
	return nothing
end


function hitmiss(l::Lininterp)
	return (l.hits,l.miss)
end


function getValue(l::Lininterp,x::Vector{Float64})
	# make sure all get evaluated
	l.ifunc = collect(1:l.nfunc)
	if l.n == 1
		eval1D(l,x)
	elseif l.n == 2
		eval2D(l,x)
	elseif l.n == 3
		eval3D(l,x)
	elseif l.n ==4
		eval4D(l,x)
	else
		warn("only up to 4D implemented so far")
	end
end
function getValue(l::Lininterp,y::Float64)
	# make sure all get evaluated
	x = [y]
	l.ifunc = collect(1:l.nfunc)
	if l.n == 1
		eval1D(l,x)
	elseif l.n == 2
		eval2D(l,x)
	elseif l.n == 3
		eval3D(l,x)
	elseif l.n ==4
		eval4D(l,x)
	else
		warn("only up to 4D implemented so far")
	end
end

# method with preallocated output vector y
# and index `which` indicating which functions to evaluate
function getValue!(y::Vector{Float64},l::Lininterp,x::Vector{Float64},which::Vector{Int})
	if maximum(which) > l.nfunc
		throw(ArgumentError("which contains a higher index than there are functions"))
	end
	if length(which) != length(y)
		throw(ArgumentError("output y must be as long as `which` index"))
	end
	# only which get evaluated
	l.ifunc = which
	if l.n == 1
		eval1D!(y,l,x)
	elseif l.n == 2
		eval2D!(y,l,x)
	elseif l.n == 3
		eval3D!(y,l,x)
	elseif l.n ==4
		eval4D!(y,l,x)
	else
		warn("only up to 4D implemented so far")
	end
end

# evaluation methods
# ==================


# finds the inf of x in it's grid for each dimension
# uses caching
# 1. remembers last evaluation, and corresponding locations in grid
# 2. finds the bracket of grid values x in is
# 3. finds position of x in the bracket
function findBracket!(l::Lininterp,x::Vector{Float64})

	# if does not have active cache, search entire interval
	if !l.hascache
		for i in 1:length(x)
			l.miss += 1
			# deal with values out of grid: set to grid bounds
			if x[i] <= l.grids[i][1]
				x[i] = l.grids[i][1]
				l.cache[i] = 1
				l.infs[i] = l.grids[i][1]
				l.sups[i] = l.grids[i][2]
				l.z[i]     = 0.0

			elseif x[i] >= l.grids[i][end]
				x[i] = l.grids[i][end]
				l.cache[i] = l.d[i]-1
				l.infs[i] = l.grids[i][end-1]
				l.sups[i] = l.grids[i][end]
				l.z[i]     = 1.0
			else
				l.cache[i] = searchsortedlast(l.grids[i],x[i],1,l.d[i]-1,Base.Forward)
				l.infs[i]  = l.grids[i][l.cache[i]]
				l.sups[i]  = l.grids[i][l.cache[i] + 1]
				l.z[i]     = (x[i] - l.infs[i]) / (l.sups[i] - l.infs[i])
			end
		end
		# now l has a cache
		l.hascache = true
	
	# else, search cache		
	else
		for i in 1:length(x)
			# deal with values out of grid: set to grid bounds
			if x[i] <= l.grids[i][1]
				l.miss += 1
				x[i] = l.grids[i][1]
				l.cache[i] = 1
				l.infs[i] = l.grids[i][1]
				l.sups[i] = l.grids[i][2]
				l.z[i]     = 0.0

			elseif x[i] >= l.grids[i][end]
				l.miss += 1
				x[i] = l.grids[i][end]
				l.cache[i] = l.d[i]-1
				l.infs[i] = l.grids[i][end-1]
				l.sups[i] = l.grids[i][end]
				l.z[i]     = 1.0
			# if x is below current lower bound of bracket, search below
			elseif x[i] < getCachedVal(l,i)
				l.miss += 1
				l.cache[i] = searchsortedlast(l.grids[i],x[i],1,l.cache[i]-1,Base.Forward)
				l.infs[i]  = l.grids[i][l.cache[i]]
				l.sups[i]  = l.grids[i][l.cache[i] + 1]
				l.z[i]     = (x[i] - l.infs[i]) / (l.sups[i] - l.infs[i])
			# if x is above current upper bound of bracket, search above
			elseif x[i] >= getNextCachedVal(l,i)
				l.miss += 1
				l.cache[i] = searchsortedlast(l.grids[i],x[i],l.cache[i]+1,l.d[i]-1,Base.Forward)
				l.infs[i]  = l.grids[i][l.cache[i]]
				l.sups[i]  = l.grids[i][l.cache[i] + 1]
				l.z[i]     = (x[i] - l.infs[i]) / (l.sups[i] - l.infs[i])
			# x is in current bracket. cool!
			else
				l.hits += 1
				l.z[i]   = (x[i] - l.infs[i]) / (l.sups[i] - l.infs[i])
			end
		end
	end
	return nothing
end


# 1D evaluation
# =============

function find1DVertices!(l::Lininterp)
	for i in l.ifunc
		l.vertex[i][1] = l.vals[i][l.cache[1]   ]	# v0
		l.vertex[i][2] = l.vals[i][l.cache[1] + 1]	# v1
	end
end

function eval1D(l::Lininterp,x::Vector{Float64})

	if length(x) != 1
		throw(ArgumentError("x needs 1 elements: one for each D in 1D!"))
	end
	if l.n != 1
		throw(ArgumentError("must supply a Lininterp with 1D"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	find1DVertices!(l)

	# build up linear combinations
	out = Float64[]
	cc = 0.0
	for i in l.ifunc
		cc = (1.0-l.z[1])* l.vertex[i][1] +   # v0
			 (    l.z[1])* l.vertex[i][2]    # v1
        push!(out,cc)
    end
    return out
end

function eval1D!(y::Vector{Float64},l::Lininterp,x::Vector{Float64})

	if length(x) != 1
		throw(ArgumentError("x needs 1 elements: one for each D in 1D!"))
	end
	if l.n != 1
		throw(ArgumentError("must supply a Lininterp with 1D"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	find1DVertices!(l)

	# build up linear combinations
	yi = 0
	for i in l.ifunc
		yi += 1
		cc = (1.0-l.z[1])* l.vertex[i][1] +   # v0
			 (    l.z[1])* l.vertex[i][2]    # v1
        y[yi] = cc
    end
    return nothing
end

# 2D evaluation
# =============

function find2DVertices!(l::Lininterp)
	for i in l.ifunc
		l.vertex[i][1] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]-1)]	# v00
		l.vertex[i][2] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1)]	# v10
		l.vertex[i][3] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]  )]	# v01
		l.vertex[i][4] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]  )]	# v11
	end
end


function eval2D(l::Lininterp,x::Array{Float64,1})

	if length(x) != 2
		throw(ArgumentError("x needs 2 elements: one for each D in 2D!"))
	end
	if l.n != 2
		throw(ArgumentError("must supply a Lininterp with 2D"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	find2DVertices!(l)

	# build up linear combinations
	out = Float64[]
	cc = 0.0
	for i in l.ifunc
		cc = (1.0-l.z[1])*(1.0-l.z[2]) * l.vertex[i][1] +   # v00
			 (    l.z[1])*(1.0-l.z[2]) * l.vertex[i][2] +   # v10
			 (1.0-l.z[1])*(    l.z[2]) * l.vertex[i][3] +   # v01
	         (    l.z[1])*(    l.z[2]) * l.vertex[i][4]     # v11
        push!(out,cc)
    end
    return out
end

function eval2D!(y::Vector{Float64},l::Lininterp,x::Array{Float64,1})

	if length(x) != 2
		throw(ArgumentError("x needs 2 elements: one for each D in 2D!"))
	end
	if l.n != 2
		throw(ArgumentError("must supply a Lininterp with 2D"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	find2DVertices!(l)

	# build up linear combinations
	yi = 0
	for i in l.ifunc
		yi += 1
		cc = (1.0-l.z[1])*(1.0-l.z[2]) * l.vertex[i][1] +   # v00
			 (    l.z[1])*(1.0-l.z[2]) * l.vertex[i][2] +   # v10
			 (1.0-l.z[1])*(    l.z[2]) * l.vertex[i][3] +   # v01
	         (    l.z[1])*(    l.z[2]) * l.vertex[i][4]     # v11
        y[yi] = cc
    end
    return nothing
end


# 3D evaluation
# =============


# finds the function values at cartesian grid of indices
# xinf, xsup, yinf, ysup, zinf, zsup
function find3DVertices!(l::Lininterp)
	for i in l.ifunc
		l.vertex[i][1] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]-1))]	# v000
		l.vertex[i][2] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]-1))]	# v100
		l.vertex[i][3] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]-1))]	# v010
		l.vertex[i][4] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]  ))]	# v001
		l.vertex[i][5] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]-1))]	# v110
		l.vertex[i][6] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]  ))]	# v101
		l.vertex[i][7] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]  ))]	# v011
		l.vertex[i][8] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]  ))]	# v111
	end
end

function eval3D(l::Lininterp,x::Array{Float64,1})

	if length(x) != 3
		throw(ArgumentError("x needs 3 elements: one for each D in 3D!"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)
	find3DVertices!(l)

	# build up linear combinations
	out = Float64[]
	cc = 0.0
	for i in l.ifunc
		cc = (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[i][1] +   # v000
			 (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[i][2] +   # v100
			 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3]) * l.vertex[i][3] +   # v010
	         (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3]) * l.vertex[i][4] +   # v001
	         (    l.z[1])*(    l.z[2])*(1.0-l.z[3]) * l.vertex[i][5] +   # v110
	         (    l.z[1])*(1.0-l.z[2])*(    l.z[3]) * l.vertex[i][6] +   # v101
	         (1.0-l.z[1])*(    l.z[2])*(    l.z[3]) * l.vertex[i][7] +   # v011
	         (    l.z[1])*(    l.z[2])*(    l.z[3]) * l.vertex[i][8]     # v111
        push!(out,cc)
	end
	return out
end

function eval3D!(y::Vector{Float64},l::Lininterp,x::Array{Float64,1})

	if length(x) != 3
		throw(ArgumentError("x needs 3 elements: one for each D in 3D!"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)
	find3DVertices!(l)

	# build up linear combinations
	yi = 0
	for i in l.ifunc
		yi += 1
		cc = (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[i][1] +   # v000
			 (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[i][2] +   # v100
			 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3]) * l.vertex[i][3] +   # v010
	         (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3]) * l.vertex[i][4] +   # v001
	         (    l.z[1])*(    l.z[2])*(1.0-l.z[3]) * l.vertex[i][5] +   # v110
	         (    l.z[1])*(1.0-l.z[2])*(    l.z[3]) * l.vertex[i][6] +   # v101
	         (1.0-l.z[1])*(    l.z[2])*(    l.z[3]) * l.vertex[i][7] +   # v011
	         (    l.z[1])*(    l.z[2])*(    l.z[3]) * l.vertex[i][8]     # v111
        y[yi] = cc
	end
	return nothing
end


# 4D evaluation
# =============

function find4DVertices!(l::Lininterp)
	for i in l.ifunc
		l.vertex[i][1]  = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v0000
		l.vertex[i][2]  = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v1000
		l.vertex[i][3]  = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v0100
		l.vertex[i][4]  = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v0010
		l.vertex[i][5]  = l.vals[i][l.cache[1] +   + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v0001
		l.vertex[i][6]  = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v1100
		l.vertex[i][7]  = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v1010
		l.vertex[i][8]  = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v1001
		l.vertex[i][9]  = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v0110
		l.vertex[i][10] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v0101
		l.vertex[i][11] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v0011
		l.vertex[i][12] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v1110
		l.vertex[i][13] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v1101
		l.vertex[i][14] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v1011
		l.vertex[i][15] = l.vals[i][l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v0111
		l.vertex[i][16] = l.vals[i][l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v1111
	end
end


function eval4D(l::Lininterp,x::Array{Float64,1})

	if length(x) != 4
		throw(ArgumentError("x needs 4 elements: one for each D in 4D!"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	find4DVertices!(l)

	wgts = hcat((1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]),
	(    l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]),
	(    l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]))

	# build up linear combinations
	m = length(l.ifunc)
	out = zeros(m)
	for i in 1:m 
		# cc = (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][1]  +   # v0000
		# 	 (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][2]  +   # v1000
		# 	 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][3]  +   # v0100
	 #         (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][4]  +   # v0010
	 #         (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][5]  +   # v0001
	 #         (    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][6]  +   # v1100
	 #         (    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][7]  +   # v1010
	 #         (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][8]  +   # v1001
		#      (1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][9]  +   # v0110
		# 	 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][10] +   # v0101
		# 	 (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][11] +   # v0011
	 #         (    l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][12] +   # v1110
	 #         (    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][13] +   # v1101
	 #         (    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][14] +   # v1011
	 #         (1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][15] +   # v0111
	 #         (    l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][16]     # v1111
        # push!(out,cc)
      	out[i] = (wgts * l.vertex[l.ifunc[i]])[1]
	end
	return out
end

function eval4D!(y::Vector{Float64},l::Lininterp,x::Array{Float64,1})

	if length(x) != 4
		throw(ArgumentError("x needs 4 elements: one for each D in 4D!"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	find4DVertices!(l)
	wgts = hcat((1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]),
	(    l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]),
	(    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]),
	(    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]),
	(1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]),
	(    l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]))

	# build up linear combinations
	yi = 0
	for i in l.ifunc
		yi += 1
		# cc = (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][1]  +   # v0000
		# 	 (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][2]  +   # v1000
		# 	 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][3]  +   # v0100
	 #         (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][4]  +   # v0010
	 #         (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][5]  +   # v0001
	 #         (    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[i][6]  +   # v1100
	 #         (    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][7]  +   # v1010
	 #         (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][8]  +   # v1001
		#      (1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][9]  +   # v0110
		# 	 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][10] +   # v0101
		# 	 (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][11] +   # v0011
	 #         (    l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[i][12] +   # v1110
	 #         (    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[i][13] +   # v1101
	 #         (    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][14] +   # v1011
	 #         (1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][15] +   # v0111
	 #         (    l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[i][16]     # v1111
        y[yi] = (wgts * l.vertex[i])[1]
	end
	return nothing
end



function show(io::IO, l::Lininterp)
	print(io,"linear interpolation object\n")
	print(io,"dimensions: $(l.d)\n")
	print(io,"approximates $(l.nfunc) functions at given point\n")
	print(io,"has active cache: $(l.hascache)\n" )
	print(io,"currently evaluates functions $(l.ifunc)\n" )
	if l.hits+l.miss == 0
		print(io,"accelerator hits 0% of attemps\n")
	else
		print(io,"accelerator hits $(round(Int,100*l.hits/(l.hits+l.miss))) % of attemps\n")
	end
end




