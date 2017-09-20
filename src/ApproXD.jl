


module ApproXD

	import Base.show, Base.convert

	using StatsBase: sample
    using PyPlot

	# load files
	include("bspline.jl")
	include("approx.jl")
	include("fspacexd.jl")
	include("lininterp.jl")
	include("plotting.jl")

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
end


