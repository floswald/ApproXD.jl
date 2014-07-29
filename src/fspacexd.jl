

type FSpaceXD
	ndim       :: Integer                                  # number of dimensions
	ndiscrete  :: Integer                                  # number of discrete choices
	coeff      :: Array{Float64}		                   # coefficient vector/matrix: ncont x ndiscrete, where ndiscrete can be 1
	idx        :: Integer				                   # idx current col index of coeff (i.e. which discrete state)
	basis      :: Dict{Integer,BSpline}	                   # dict with spline objects
	curr_basis :: Dict{Int,SparseMatrixCSC{Float64,Int64}} # dict with current spline evals

	function FSpaceXD(ndim::Integer,coeff::Array{Float64},bs::Dict{Integer,BSpline})
		if ndim != length(bs)
			throw(ArgumentError("number of basis functions must be equal to ndim"))
		end
		cols = 1
		for (k,v) in bs
			cols *= getNumCoefs(bs[k])
		end
		if cols != size(coeff,1)
			throw(ArgumentError("Dimension mismatch: size(coeff,1) must be equal to the product of all spline coefficients"))
		end
		cb = Dict{Int,SparseMatrixCSC{Float64,Int64}}()
		for i in ndim
			cb[i] = spzeros(getNumCoefs(bs[i]),1)
		end

		# number of dchoices
		if ndims(coeff) > 1
			ndiscrete = size(coeff,2)
			idx = 1
		else
			ndiscrete = 1
			idx = 0
		end
		return new(ndim,ndiscrete,coeff,idx,bs,cb)
	end
end

function setindex!(fx::FSpaceXD,i::Integer)
	fx.idx = i
end

function getValue(x::Array{Float64,1},fx::FSpaceXD)

	# 0) check length of x
	if length(x) != fx.ndim
		throw(ArgumentError("length of x must be equal to fx.ndim"))
	end

	# 1) evaluate each spline[j] at x[j]
	for i in 1:fx.ndim
		fx.curr_basis[i] = getBasis(x[i],fx.basis[i])
	end

	# 2) use evalTensorXX to evaluate at correct coeff vector (in case there are many)
	if fx.ndiscrete > 1
		evalTensor(fx.curr_basis,fx.coeff[:,fx.idx])
	else
		evalTensor(fx.curr_basis,fx.coeff)
	end
end
