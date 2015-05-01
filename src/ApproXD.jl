


module ApproXD

	import Base.show, Base.convert

	# load files
	include("bspline.jl")
	include("approx.jl")
	include("fspacexd.jl")
	include("lininterp.jl")

	if Sys.OS_NAME == :Darwin
	   using PyPlot
	   include("plotting.jl")
	end

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


