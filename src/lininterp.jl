

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
		if n > 4
			throw(ArgumentError("currently only up to 4D implemented"))
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


function getValue(l::lininterp,x::Vector{Float64})
	if l.n == 2
		eval2D(l,x)
	elseif l.n == 3
		eval3D(l,x)
	elseif l.n ==4
		eval4D(l,x)
	else
		warn("only up to 4D implemented so far")
	end
end

function eval2D(l::lininterp,x::Array{Float64,1})

	if length(x) != 2
		throw(ArgumentError("x needs 2 elements: one for each D in 2D!"))
	end
	if l.n != 2
		throw(ArgumentError("must supply a lininterp with 2D"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	if l.hitnow
		# use the same function values as in the last iteration!
		# therefore do nothing
	else
		# if missed, need to find the function values on the new bracket
		find2DVertices!(l)
	end

	# build up linear combinations
	cc = (1.0-l.z[1])*(1.0-l.z[2]) * l.vertex[1] +   # v00
		 (    l.z[1])*(1.0-l.z[2]) * l.vertex[2] +   # v10
		 (1.0-l.z[1])*(    l.z[2]) * l.vertex[3] +   # v01
         (    l.z[1])*(    l.z[2]) * l.vertex[4]     # v11
	return cc
end

function eval3D(l::lininterp,x::Array{Float64,1})

	if length(x) != 3
		throw(ArgumentError("x needs 3 elements: one for each D in 3D!"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	if l.hitnow
		# use the same function values as in the last iteration!
		# therefore do nothing
	else
		# if missed, need to find the function values on the new bracket
		find3DVertices!(l)
	end

	# build up linear combinations
	cc = (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3]) * l.vertex[1] +   # v000
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

function eval4D(l::lininterp,x::Array{Float64,1})

	if length(x) != 4
		throw(ArgumentError("x needs 4 elements: one for each D in 4D!"))
	end

	# find in which bracket of grid values x is in.
	# using cached values
	# finds the relative position of x inside the bracket
	findBracket!(l,x)

	# get the function values on the box
	if l.hitnow
		# use the same function values as in the last iteration!
		# therefore do nothing
	else
		# if missed, need to find the function values on the new bracket
		find4DVertices!(l)
	end

	# build up linear combinations
	cc = (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[1]  +   # v0000
		 (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[2]  +   # v1000
		 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[3]  +   # v0100
         (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[4]  +   # v0010
         (1.0-l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[5]  +   # v0001
         (    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(1.0-l.z[4]) * l.vertex[6]  +   # v1100
         (    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[7]  +   # v1010
         (    l.z[1])*(1.0-l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[8]  +   # v1001
	     (1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[9]  +   # v0110
		 (1.0-l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[10] +   # v0101
		 (1.0-l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[11] +   # v0011
         (    l.z[1])*(    l.z[2])*(    l.z[3])*(1.0-l.z[4]) * l.vertex[12] +   # v1110
         (    l.z[1])*(    l.z[2])*(1.0-l.z[3])*(    l.z[4]) * l.vertex[13] +   # v1101
         (    l.z[1])*(1.0-l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[14] +   # v1011
         (1.0-l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[15] +   # v0111
         (    l.z[1])*(    l.z[2])*(    l.z[3])*(    l.z[4]) * l.vertex[16]     # v1111

	return cc
end

function find2DVertices!(l::lininterp)
	l.vertex[1] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]-1)]	# v00
	l.vertex[2] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1)]	# v10
	l.vertex[3] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]  )]	# v01
	l.vertex[4] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]  )]	# v11
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

function find4DVertices!(l::lininterp)
	l.vertex[1]  = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v0000
	l.vertex[2]  = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v1000
	l.vertex[3]  = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v0100
	l.vertex[4]  = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v0010
	l.vertex[5]  = l.vals[l.cache[1] +   + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v0001
	l.vertex[6]  = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]-1)))]	# v1100
	l.vertex[7]  = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v1010
	l.vertex[8]  = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v1001
	l.vertex[9]  = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v0110
	l.vertex[10] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v0101
	l.vertex[11] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v0011
	l.vertex[12] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]-1)))]	# v1110
	l.vertex[13] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]-1) + l.d[3]*(l.cache[4]  )))]	# v1101
	l.vertex[14] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]-1 + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v1011
	l.vertex[15] = l.vals[l.cache[1] +     l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v0111
	l.vertex[16] = l.vals[l.cache[1] + 1 + l.d[1]*(l.cache[2]   + l.d[2]*((l.cache[3]  ) + l.d[3]*(l.cache[4]  )))]	# v1111
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
		hits = 0
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
				hits += 1
				l.z[i]   = (x[i] - l.infs[i]) / (l.sups[i] - l.infs[i])
			end
		end
		l.hitnow = hits == length(x)
	end
	return nothing
end

function show(io::IO, l::lininterp)
	print(io,"linear interpolation object with\n")
	print(io,"dimensions: $(l.d)\n")
	print(io,"has active cache: $(l.hascache)\n" )
	if l.hits+l.miss == 0
		print(io,"accelerator hits 0% of attemps\n")
	else
		print(io,"accelerator hits $(round(100*l.hits/(l.hits+l.miss),2)) % of attemps\n")
	end
end




