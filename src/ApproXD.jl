module ApproXD

import Base.show, Base.convert

	using StatsBase: sample
	using SparseArrays
	using LinearAlgebra


	# load files
	include("bspline.jl")
	include("approx.jl")
	include("fspacexd.jl")
	include("lininterp.jl")
	# couldn't get PyPlot to install properly on appveyor.
	# so windows users can't do the plots. there not that important anyway.
	# include("plotting.jl")

	export BSpline,
		   Lininterp,
	       show,
	       getBasis,
	       getNumCoefs,
	       getNumKnots,
	       getTensorCoef,
	       FSpaceXD,
	       getValue,
	       setValue,
	       setindex!,
	       getValue!,
	       getValue,
	       hitmiss,
	       resetCache!,
	       setGrid!

end # module
