





# main programme
home = ENV["HOME"]
cd("$home/git/BSplines.jl")

include("src/BSplines.jl")
include("test/test_lininterp.jl")

# run individual tests
include("test/test_basics.jl")
include("test/test_approx.jl")
include("test/test_FSpaceXD.jl")


# run all tests: exits 
include("test/runtests.jl")

