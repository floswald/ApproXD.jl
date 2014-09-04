

module test_lininterp

using BSplines, FactCheck

facts("constructor for lininterp") do

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, linspace(1.0,3,3))
	push!(gs, linspace(2.0,3,4))
	push!(gs, linspace(-1,3.0,5))

	l = BSplines.lininterp(vs,gs)

	@fact BSplines.getDims(l) => [3,4,5]


	# errors


end






end
