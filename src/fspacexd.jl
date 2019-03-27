

mutable struct FSpaceXD
	ndim       :: Integer                                  # number of dimensions
	coeff      :: Array{Float64}		                   # coefficient vector/
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

		return new(ndim,coeff,bs,cb)
	end
end

function getValue(fx::FSpaceXD,x::Array{Float64,1})

	# 0) check length of x
	if length(x) != fx.ndim
		throw(ArgumentError("length of x must be equal to fx.ndim"))
	end

	# 1) evaluate each spline[j] at x[j]
	for i in 1:fx.ndim
		fx.curr_basis[i] = getBasis(x[i],fx.basis[i])
	end

	# 2) use evalTensorXX to evaluate at correct coeff vector (in case there are many)
	evalTensor(fx.curr_basis,fx.coeff)
end


function show(io::IO,f::FSpaceXD)
	print(io,"FSpaceXD with $(f.ndim) dimensions\n")
end

