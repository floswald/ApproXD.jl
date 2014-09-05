

# multidimensional linear interpolator
# with cache accelerator
type lininterp

	d        :: Array{Int}    # number of points in each dimension
	n        :: Int 	        # number  of dims
	hascache :: Bool 	        # has active cache
	cache    :: Array{Int}	# current index of lower bound of bracket containing x
	infs     :: Array{Float64}	    # current value of inf: lower bound of bracket in xgrid
	sups     :: Array{Float64}	   	# current value of sup: upper bound of bracket in xgrid
	z        :: Array{Float64}    	# current position in bracket: z = (x-x[inf]) / (x[sup] - x[inf])
	hits     :: Int 			# count of cache hits
	miss     :: Int 			# count of cache misses
	hitnow   :: Bool 			# whether current eval is a hit
	grids    :: Array{Array{Float64}} 	# grids for each dimension
	vals     :: Array{Float64} 	# array with function values on the grid
	vertex   :: Array{Float64}  # array with current vertices

	function lininterp(v::Array{Float64},g::Array{Array{Float64,1}})
		d = [size(v)...]
		n = length(d)
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
		if n > 3
			throw(ArgumentError("currently only up to 3D implemented"))
		end
		cache = ones(Int,n)
		return new(d,n,false,cache,zeros(n),zeros(n),zeros(n),0,0,false,g,v,zeros(2^n))
	end

end


getDims(l::lininterp) = l.d 
getGrids(l::lininterp) =  l.grids 
function resetCache!(l::lininterp) 
	fill!(l.cache,0)
	l.hits = 0
	l.miss = 0
	l.hascache=false
	return nothing
end
function getCache(l::lininterp)
	l.cache
end
function getCache(l::lininterp,i::Int)
	l.cache[i]
end
function getCachedVal(l::lininterp,i::Int)
	l.grids[i][l.cache[i]]
end
function getNextCachedVal(l::lininterp,i::Int)
	l.grids[i][l.cache[i] + 1]
end

function setVals(l::lininterp,v::Array)
	if size(v) != d
		throw(ArgumentError("size(v) must be equal to size d"))
	end
	vals = v
	return nothing
end

function hitmiss(l::lininterp)
	return (l.hits,l.miss)
end


function eval(l::lininterp,x::Vector{Float64})
	if l.n == 3
		eval3D(l,x)
	else
		warn("only 3D implemented so far")
	end
end

function eval3D(l::lininterp,z::Array{Float64,1})

	if length(z) != 3
		throw(ArgumentError("x needs 3 elements"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,z)

	# get the function values on the box
	# TODO why does this not work?
	# if l.hitnow
		# use the same function values as in the last iteration!
		# therefore do nothing
	# else
		# if missed, need to find the function values on the new bracket
		find3DVertices!(l)
	# end

	# build up linear combinations
	cc =      (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[1] +   # v000
		      (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[2] +   # v100
		      (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3]) * l.vertex[3] +   # v010
              (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3]) * l.vertex[4] +   # v001
              (    l.z[1])*(    l.z[2])*(1.0-l.z[3]) * l.vertex[5] +   # v110
              (    l.z[1])*(1.0-l.z[2])*(    l.z[3]) * l.vertex[6] +   # v101
              (1.0-l.z[1])*(    l.z[2])*(    l.z[3]) * l.vertex[7] +   # v011
              (    l.z[1])*(    l.z[2])*(    l.z[3]) * l.vertex[8]     # v111


	# v = cc * l.vertex
	return cc
end

# finds the function values at cartesian grid of indices
# xinf, xsup, yinf, ysup, zinf, zsup
function find3DVertices!(l::lininterp)
	l.vertex[1] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]-1))]	# v000
	l.vertex[2] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]-1))]	# v100
	l.vertex[3] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]-1))]	# v010
	l.vertex[4] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]  ))]	# v001
	l.vertex[5] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]-1))]	# v110
	l.vertex[6] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*(l.cache[3]  ))]	# v101
	l.vertex[7] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]  ))]	# v011
	l.vertex[8] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*(l.cache[3]  ))]	# v111
end

# finds the inf of x in it's grid for each dimension
# uses caching
# 1. remembers last evaluation, and corresponding locations in grid
# 2. sets boolean hitnow if current evaluation is equal. if true, 
# can use current values in vertex (you were in exactly that bracket last time)
# 3. finds the bracket of grid values x in is
# 4. finds position of x in the bracket
function findBracket!(l::lininterp,x::Vector{Float64})

	l.hitnow = false

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
				l.hitnow = true
			end
		end
	end
	return nothing
end

function show(io::IO, l::lininterp)
	print(io,"linear interpolation object with\n")
	print(io,"number of dims: $(l.n)\n")
	print(io,"has active cache: $(l.hascache)\n" )
	print(io,"current vertices: $(l.vertex)\n" )
	if l.hits+l.miss == 0
		print(io,"accelerator hits 0% of attemps\n")
	else
		print(io,"accelerator hits $(round(100*l.hits/(l.hits+l.miss),2)) % of attemps\n")
	end
end




