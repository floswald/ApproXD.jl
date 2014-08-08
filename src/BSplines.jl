


module BSplines

	import Base.show
	using PyPlot, ArrayViews

	# load files
	include("bspline.jl")
	include("approx.jl")
	include("fspacexd.jl")

	export BSpline,
	       show,
	       getBasis,
	       getNumCoefs,
	       getNumKnots,
	       getTensorCoef,
	       FSpaceXD,
	       getValue,
	       setindex!

end


