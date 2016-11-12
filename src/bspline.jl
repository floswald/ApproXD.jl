




type BSpline
	degree   :: Int
	numKnots :: Int
	lower    :: Float64
	upper    :: Float64
	knots    :: Array{Float64,1}

	# constructor with equally spaced interior knots
	function BSpline(nKnots,deg,lb,ub)
		knots = zeros(nKnots+2*deg)
		if (nKnots < 2*(deg+1)-1)
			throw(ArgumentError("need at least 2*(deg+1) - 1 = $(2*(deg+1) - 1) knots"))
		end
		if deg < 0
			throw(ArgumentError("degree must be non-negative integer"))
		end
		if lb > ub
			throw(ArgumentError("require lb < ub"))
		end

		h = (ub - lb) / (nKnots - 1)
		for i=1:length(knots)
			if i < deg + 1
				knots[i] = lb
			elseif i > nKnots + deg 
				knots[i] = ub
			else
				knots[i] = lb + (i-deg-1)*h
			end
		end
		# sp = sparsevec(zeros(num_coefs))
		new(deg,nKnots,lb,ub,knots)
	end

	# constructor with user supplied interior knots
	function BSpline(knots::Array,deg::Integer)
		if !issorted(knots)
			throw(ArgumentError("knots must be sorted"))
		end
		lb = knots[1]
		ub = knots[end]
		# extend knot vector by degree on both sides
		knots_new = [ [lb for i=1:deg]; knots; [ub for i=1:deg] ]
		numKnots = length(knots)	# number of INTERIOR knots 
		new(deg,numKnots,lb,ub,knots_new)
	end

end


function show(io::IO, b::BSpline)
	print(io,"BSpline object with\n")
	print(io,"degree: $(b.degree)\n")
	print(io,"number of knots: $(b.numKnots)\n")
	print(io,"[lower,upper]: [$(b.lower),$(b.upper)]\n")
	print(io,"knot vector: $(b.knots)\n")
end

function getNumKnots(b::BSpline) return b.numKnots end
function getNumCoefs(b::BSpline) return b.numKnots + b.degree -1 end
function getCoefs(b::BSpline) return b.knots[ (b.degree+1):(length(b.knots)-b.degree)] end
function getDegree(b::BSpline) return b.degree end


# """
# 	compute value of the ``d + 1`` B-splines at ``x``

# returns 

# ```math
# b = \right(B_{\mu-d,d}(x),\dots,B_{\mu,d}(x)\left)^T
# ```
# """
# function getBasis2(x::Float64,b::BSpline)

# 	num_nodes = getNumCoefs(b)
# 	deg       = b.degree

# 	B = zeros(num_nodes)  # Bspline

# end
 
#
function getBasis(x::Float64,b::BSpline)
	
	num_nodes = getNumCoefs(b)
	deg       = b.degree
	
	# tmp
	d = 0.0
	e = 0.0

	# create a basis function
	# return a colvec because currntly only CSC format 
	# a 1-row sparse matrix is dense
	# bs = spzeros(num_nodes,1)
	bs = spzeros(num_nodes,1)

	# check x
	if x < b.lower
		warn("x < lb: x=$x, lb=$(b.lower). be careful!")
	elseif x > b.upper
		warn("x > ub: x=$x, ub=$(b.upper). be careful!")
	end

	# get mu s.t. knot_mu < knot_{mu+1} and x in [knot_mu, knot_{mu+1})
	# i.e. get the index of the lower knot in the active knot span

	# fix bound behaviour
	if x <= b.lower
		mu = deg+1
	elseif x>=b.upper
		mu = num_nodes
	else
		mu = searchsortedlast(b.knots,x) 
	end

	# set 0-degree basis function
	# 0-deg basis is an indicator function
	# that is 1.0 in the active knot span and 0.0 else.
	bs[mu] = 1.0

	# loop over degrees
	for k=1:deg

		# loop over basis functions
		for j in mu-k:mu

			# take care of "division by zero" issue
			# dividing by zero must return 0.0
			if j+k <= deg +1
				d = 0.0
			elseif b.knots[j+k]-b.knots[j] == 0
				d = 0.0
			else
				d = bs[j] * (x - b.knots[j]) / (b.knots[j+k]-b.knots[j])
			end

			if j+1 >= num_nodes+1
				e = 0.0
			elseif b.knots[j+k+1]-b.knots[j+1] == 0
				e = 0.0
			else
				e = bs[j+1] * (b.knots[j+k+1] - x)/(b.knots[j+k+1]-b.knots[j+1])
			end

			bs[j] = d + e

		end

	end
	return bs
end

# vector of points
function getBasis(x::Vector{Float64},b::BSpline)

	n = length(x)
	num_nodes = getNumCoefs(b)
	deg       = b.degree


	# if you evaluate a degree 1 basis at its interior knots ("coefs"), you get identity
	if (deg == 1) && (n == num_nodes) && (sum(abs(x .- getCoefs(b))) < 1e-8)
		return speye(n)
	else
	# sort x?
		
		# tmp
		d = 0.0
		e = 0.0

		# create a basis function
		# return a colvec because currntly only CSC format 
		# a 1-row sparse matrix is dense
		bs = spzeros(n,num_nodes)

		for xi in 1:n

			# check x
			if x[xi] < b.lower
				warn("x < lb: x=$(x[xi]), lb=$(b.lower). be careful!")
			elseif x[xi] > b.upper
				warn("x > ub: x=$(x[xi]), ub=$(b.upper). be careful!")
			end

			# get mu s.t. knot_mu < knot_{mu+1} and x in [knot_mu, knot_{mu+1})
			# i.e. get the index of the lower knot in the active knot span

			# fix bound behaviour
			if x[xi] <= b.lower
				mu = deg+1
			elseif x[xi]>=b.upper
				mu = num_nodes
			else
				mu = searchsortedlast(b.knots,x[xi]) 
			end

			# set 0-degree basis function
			# 0-deg basis is an indicator function
			# that is 1.0 in the active knot span and 0.0 else.
			bs[xi,mu] = 1.0

			# loop over degrees
			for k=1:deg

				# loop over basis functions
				for j in mu-k:mu

					# take care of "division by zero" issue
					# dividing by zero must return 0.0
					if j+k <= deg +1
						d = 0.0
					elseif b.knots[j+k]-b.knots[j] == 0
						d = 0.0
					else
						d = bs[xi,j] * (x[xi] - b.knots[j]) / (b.knots[j+k]-b.knots[j])
					end

					if j+1 >= num_nodes+1
						e = 0.0
					elseif b.knots[j+k+1]-b.knots[j+1] == 0
						e = 0.0
					else
						e = bs[xi,j+1] * (b.knots[j+k+1] - x[xi])/(b.knots[j+k+1]-b.knots[j+1])
					end
					bs[xi,j] = d + e
				end

			end
		end
		return bs

	end
	
end