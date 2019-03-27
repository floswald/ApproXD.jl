
module test_lin
using ApproXD
using Test
# using SparseArrays
# using LinearAlgebra

@testset "constructor for 1D Lininterp on a single function" begin

	vs = rand(4)
	gs = Array{Float64,1}[]
	push!(gs, range(2.0, stop = 3, length = 4))

	l = Lininterp(vs,gs)

	@test ApproXD.getDims(l) == [4]
	@test l.n == 1
	@test l.nfunc == 1
	@test l.ifunc == [1]
	@test l.infs == zeros(1)
	@test l.sups == zeros(1)
	@test l.hits == 0
	@test l.miss == 0
	@test l.vals[1] == vs
	@test l.grids == gs
	@test ApproXD.getCache(l) == [1]

	# errors
	pop!(gs)   # missing grid
	@test_throws ArgumentError Lininterp(vs,gs)   
	push!(gs, range(-1,stop =3.0, length=5)) 	# wrong length
	@test_throws ArgumentError Lininterp(vs,gs)
	pop!(gs)
	push!(gs, [1.1,0,2,3]) 	# not sorted
	@test_throws ArgumentError Lininterp(vs,gs)

end

@testset "constructor for 2D Lininterp on a single function" begin

	vs = rand(3,4)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0, stop = 3, length = 4))

	l = Lininterp(vs,gs)

	@test ApproXD.getDims(l) == [3,4]
	@test l.n == 2
	@test l.nfunc == 1
	@test l.ifunc == [1]
	@test l.infs == zeros(2)
	@test l.sups == zeros(2)
	@test l.hits == 0
	@test l.miss == 0
	@test l.vals[1] == vs
	@test l.grids == gs
	@test ApproXD.getCache(l) == [1,1]

	# errors
	pop!(gs)   # missing grid
	@test_throws ArgumentError Lininterp(vs,gs)   
	push!(gs, range(-1,stop =3.0, length=5)) 	# wrong length
	@test_throws ArgumentError Lininterp(vs,gs)
	pop!(gs)
	push!(gs, [1.1,0,2,3]) 	# not sorted
	@test_throws ArgumentError Lininterp(vs,gs)

end

@testset "constructor for 2D Lininterp on two functions" begin

	vs1 = rand(3,4)
	vs2 = rand(3,4)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0, stop = 3, length = 4))

	l = Lininterp(vs1,vs2,gs)

	@test ApproXD.getDims(l) == [3,4]
	@test l.n == 2
	@test l.nfunc == 2
	@test l.ifunc == [1,2]
	@test l.infs == zeros(2)
	@test l.sups == zeros(2)
	@test l.hits == 0
	@test l.miss == 0
	@test l.vals[1] == vs1
	@test l.vals[2] == vs2
	@test l.grids == gs
	@test ApproXD.getCache(l) == [1,1]

	# errors
	pop!(gs)   # missing grid
	@test_throws UndefVarError Lininterp(vs,gs)   
	push!(gs, range(-1,stop =3.0, length=5)) 	# wrong length
	@test_throws UndefVarError Lininterp(vs,gs)
	pop!(gs)
	push!(gs, [1.1,0,2,3]) 	# not sorted
	@test_throws UndefVarError Lininterp(vs,gs)

end

@testset "constructor for 3D Lininterp" begin

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0, stop = 3, length = 4))
	push!(gs, range(-1,stop =3.0, length=5))

	l = Lininterp(vs,gs)

	@test ApproXD.getDims(l) == [3,4,5]
	@test l.n == 3
	@test l.ifunc == [1]
	@test l.infs == zeros(3)
	@test l.sups == zeros(3)
	@test l.hits == 0
	@test l.miss == 0
	@test l.vals[1] == vs
	@test l.grids == gs
	@test ApproXD.getCache(l) == [1,1,1]

	# errors
	pop!(gs)   # missing grid
	@test_throws ArgumentError Lininterp(vs,gs)   
	push!(gs, range(-1, stop = 3.0, length = 4)) 	# wrong length
	@test_throws ArgumentError Lininterp(vs,gs)
	pop!(gs)
	push!(gs, [1.1,0,2,3,4]) 	# not sorted
	@test_throws ArgumentError Lininterp(vs,gs)

end

@testset "constructor for 4D Lininterp" begin

	vs = rand(3,4,5,7)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0, stop = 3, length = 4))
	push!(gs, range(-1,stop =3.0, length=5))
	push!(gs, range(-0.1,stop = 4.8,length = 7))

	l = Lininterp(vs,gs)

	@test ApproXD.getDims(l) == [3,4,5,7]
	@test l.n == 4
	@test l.infs == zeros(4)
	@test l.sups == zeros(4)
	@test l.hits == 0
	@test l.miss == 0
	@test l.vals[1] == vs
	@test l.grids == gs
	@test ApproXD.getCache(l) == [1,1,1,1]

	# errors
	pop!(gs)   # missing grid
	@test_throws ArgumentError Lininterp(vs,gs)   
	push!(gs, range(-1, stop = 3.0, length = 4)) 	# wrong length
	@test_throws ArgumentError Lininterp(vs,gs)
	pop!(gs)
	push!(gs, [1.1,0,2,3,4,5,6]) 	# not sorted
	@test_throws ArgumentError Lininterp(vs,gs)

end


@testset "testing Cache functions" begin

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0, stop = 3, length = 4))
	push!(gs, range(-1,stop =3.0, length=5))

	l = Lininterp(vs,gs)
	myc = [2,1,3]
	l.cache = myc	# manually set a cache

	for i in 1:length(gs)
		@test ApproXD.getCache(l,i) == myc[i]
		@test ApproXD.getCachedVal(l,i) == gs[i][myc[i]]
		@test ApproXD.getNextCachedVal(l,i) == gs[i][myc[i]+1]
	end


end

@testset "testing findBracket!(x) for 1D" begin

	vs = rand(5)
	gs = Array{Float64,1}[]
	push!(gs, range(-1,stop =3.0, length=5))

	l = Lininterp(vs,gs)
	@test l.hascache == false

	x = [2.8] 	# == brackets 4
	ApproXD.findBracket!(l,x)

	# brackets is index of lower bound of each bracket
	brackets = [4]

	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (0,1)
	@test l.infs[1] == gs[1][brackets[1]]
	@test l.sups[1] == gs[1][brackets[1]+1]

	@test l.hascache == true
	ApproXD.findBracket!(l,x)
	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (1,1)

	# value x1 and x2 out of bounds
	x = [3.8] 	# == brackets 4
	ApproXD.findBracket!(l,x)

	brackets = [4]
	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (1,2)
	@test l.infs[1] == gs[1][end-1]
	@test l.sups[1] == gs[1][end]

	# side-effect: x is forced into bounds!
	# may or may not be desirable.
	@test x[1] == l.sups[1]
end



@testset "testing findBracket!(x) for 3D" begin

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0,stop = 6.0,length = 4))
	push!(gs, range(-1,stop =3.0, length=5))

	l = Lininterp(vs,gs)
	@test l.hascache == false

	x = [1.1,3.5,2.8] 	# == brackets 1,2,4
	ApproXD.findBracket!(l,x)

	# brackets is index of lower bound of each bracket
	brackets = [1,2,4]

	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (0,3)
	for i in 1:3
		@test l.infs[i] == gs[i][brackets[i]]
		@test l.sups[i] == gs[i][brackets[i]+1]
	end

	@test l.hascache == true
	ApproXD.findBracket!(l,x)
	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (3,3)

	# value x1 and x2 out of bounds
	x = [0.9,7.5,2.8] 	# == brackets 1,3,4
	ApproXD.findBracket!(l,x)

	brackets = [1,3,4]
	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (4,5)
	@test l.infs[1] == gs[1][1]
	@test l.infs[2] == gs[2][end-1]
	@test l.infs[3] == gs[3][brackets[3]]
	@test l.sups[1] == gs[1][2]
	@test l.sups[2] == gs[2][end]
	@test l.sups[3] == gs[3][brackets[3]+1]

	# side-effect: x is forced into bounds!
	# may or may not be desirable.
	@test x[1] == l.infs[1]
	@test x[2] == l.sups[2]
end

@testset "testing findBracket!(x) for 4D" begin

	vs = rand(3,4,5,7)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0,stop = 6.0,length = 4))
	push!(gs, range(-1,stop =3.0, length=5))
	push!(gs, range(-2,stop = 5.0,length = 7))

	l = Lininterp(vs,gs)
	@test l.hascache == false

	x = [1.1,3.5,2.8,3.9] 	# == brackets 1,2,4,6
	ApproXD.findBracket!(l,x)

	# brackets is index of lower bound of each bracket
	brackets = [1,2,4,6]

	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (0,4)
	for i in 1:4
		@test l.infs[i] == gs[i][brackets[i]]
		@test l.sups[i] == gs[i][brackets[i]+1]
	end

	@test l.hascache == true
	ApproXD.findBracket!(l,x)
	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (4,4)

	# value x1 and x2 out of bounds
	x = [0.9,7.5,2.8,3.9] 	# == brackets 1,3,4,6
	ApproXD.findBracket!(l,x)

	brackets = [1,3,4,6]
	@test ApproXD.getCache(l) == brackets
	@test ApproXD.hitmiss(l) == (6,6)
	@test l.infs[1] == gs[1][1]
	@test l.infs[2] == gs[2][end-1]
	@test l.infs[3] == gs[3][brackets[3]]
	@test l.infs[4] == gs[4][brackets[4]]
	@test l.sups[1] == gs[1][2]
	@test l.sups[2] == gs[2][end]
	@test l.sups[3] == gs[3][brackets[3]+1]
	@test l.sups[4] == gs[4][brackets[4]+1]

	# side-effect: x is forced into bounds!
	# may or may not be desirable.
	@test x[1] == l.infs[1]
	@test x[2] == l.sups[2]
end





@testset "testing find3DVertices!(l) for 1 function" begin
	
	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0,stop = 6.0,length = 4))
	push!(gs, range(-1,stop =3.0, length=5))

	l = Lininterp(vs,gs)
	x = [1.1,3.5,2.8] 	# == brackets 1,2,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,3)
	ApproXD.find3DVertices!(l)

	@test length(l.vertex[1]) == 8
	@test l.vertex[1][1] == vs[1,2,4]  # v000   
	@test l.vertex[1][2] == vs[2,2,4]  # v100
	@test l.vertex[1][3] == vs[1,3,4]  # v010
	@test l.vertex[1][4] == vs[1,2,5]  # v001
	@test l.vertex[1][5] == vs[2,3,4]  # v110
	@test l.vertex[1][6] == vs[2,2,5]  # v101
	@test l.vertex[1][7] == vs[1,3,5]  # v011
	@test l.vertex[1][8] == vs[2,3,5]  # v111

	x = [3,6,3.0] 	# == brackets 2,3,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,6)
	ApproXD.find3DVertices!(l)

	@test l.vertex[1][1] == vs[2,3,4]  # v000   
	@test l.vertex[1][2] == vs[3,3,4]  # v100
	@test l.vertex[1][3] == vs[2,4,4]  # v010
	@test l.vertex[1][4] == vs[2,3,5]  # v001
	@test l.vertex[1][5] == vs[3,4,4]  # v110
	@test l.vertex[1][6] == vs[3,3,5]  # v101
	@test l.vertex[1][7] == vs[2,4,5]  # v011
	@test l.vertex[1][8] == vs[3,4,5]  # v111

	# out of bounds: no hits by convention
	x = [-1,6.5,-3.0] 	# == brackets 1,3,1
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,9)
	ApproXD.find3DVertices!(l)

	@test l.vertex[1][1] == vs[1,3,1]  # v000   
	@test l.vertex[1][2] == vs[2,3,1]  # v100
	@test l.vertex[1][3] == vs[1,4,1]  # v010
	@test l.vertex[1][4] == vs[1,3,2]  # v001
	@test l.vertex[1][5] == vs[2,4,1]  # v110
	@test l.vertex[1][6] == vs[2,3,2]  # v101
	@test l.vertex[1][7] == vs[1,4,2]  # v011
	@test l.vertex[1][8] == vs[2,4,2]  # v111
end


@testset "testing find3DVertices!(l) for 2 functions" begin
	
	vs1 = rand(3,4,5)
	vs2 = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0,stop = 6.0,length = 4))
	push!(gs, range(-1,stop =3.0, length=5))

	l = Lininterp(vs1,vs2,gs)
	x = [1.1,3.5,2.8] 	# == brackets 1,2,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,3)
	ApproXD.find3DVertices!(l)

	@test length(l.vertex[1]) == 8
	@test l.vertex[1][1] == vs1[1,2,4]  # v000   
	@test l.vertex[1][2] == vs1[2,2,4]  # v100
	@test l.vertex[1][3] == vs1[1,3,4]  # v010
	@test l.vertex[1][4] == vs1[1,2,5]  # v001
	@test l.vertex[1][5] == vs1[2,3,4]  # v110
	@test l.vertex[1][6] == vs1[2,2,5]  # v101
	@test l.vertex[1][7] == vs1[1,3,5]  # v011
	@test l.vertex[1][8] == vs1[2,3,5]  # v111

	@test l.vertex[2][1] == vs2[1,2,4]  # v000   
	@test l.vertex[2][2] == vs2[2,2,4]  # v100
	@test l.vertex[2][3] == vs2[1,3,4]  # v010
	@test l.vertex[2][4] == vs2[1,2,5]  # v001
	@test l.vertex[2][5] == vs2[2,3,4]  # v110
	@test l.vertex[2][6] == vs2[2,2,5]  # v101
	@test l.vertex[2][7] == vs2[1,3,5]  # v011
	@test l.vertex[2][8] == vs2[2,3,5]  # v111

	x = [3,6,3.0] 	# == brackets 2,3,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,6)
	ApproXD.find3DVertices!(l)

	@test l.vertex[1][1] == vs1[2,3,4]  # v000   
	@test l.vertex[1][2] == vs1[3,3,4]  # v100
	@test l.vertex[1][3] == vs1[2,4,4]  # v010
	@test l.vertex[1][4] == vs1[2,3,5]  # v001
	@test l.vertex[1][5] == vs1[3,4,4]  # v110
	@test l.vertex[1][6] == vs1[3,3,5]  # v101
	@test l.vertex[1][7] == vs1[2,4,5]  # v011
	@test l.vertex[1][8] == vs1[3,4,5]  # v111

	@test l.vertex[2][1] == vs2[2,3,4]  # v000   
	@test l.vertex[2][2] == vs2[3,3,4]  # v100
	@test l.vertex[2][3] == vs2[2,4,4]  # v010
	@test l.vertex[2][4] == vs2[2,3,5]  # v001
	@test l.vertex[2][5] == vs2[3,4,4]  # v110
	@test l.vertex[2][6] == vs2[3,3,5]  # v101
	@test l.vertex[2][7] == vs2[2,4,5]  # v011
	@test l.vertex[2][8] == vs2[3,4,5]  # v111

	# out of bounds: no hits by convention
	x = [-1,6.5,-3.0] 	# == brackets 1,3,1
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,9)
	ApproXD.find3DVertices!(l)

	@test l.vertex[1][1] == vs1[1,3,1]  # v000   
	@test l.vertex[1][2] == vs1[2,3,1]  # v100
	@test l.vertex[1][3] == vs1[1,4,1]  # v010
	@test l.vertex[1][4] == vs1[1,3,2]  # v001
	@test l.vertex[1][5] == vs1[2,4,1]  # v110
	@test l.vertex[1][6] == vs1[2,3,2]  # v101
	@test l.vertex[1][7] == vs1[1,4,2]  # v011
	@test l.vertex[1][8] == vs1[2,4,2]  # v111

	@test l.vertex[2][1] == vs2[1,3,1]  # v000   
	@test l.vertex[2][2] == vs2[2,3,1]  # v100
	@test l.vertex[2][3] == vs2[1,4,1]  # v010
	@test l.vertex[2][4] == vs2[1,3,2]  # v001
	@test l.vertex[2][5] == vs2[2,4,1]  # v110
	@test l.vertex[2][6] == vs2[2,3,2]  # v101
	@test l.vertex[2][7] == vs2[1,4,2]  # v011
	@test l.vertex[2][8] == vs2[2,4,2]  # v111
end




@testset "testing find4DVertices!(l) single function" begin
	
	vs = rand(3,4,5,7)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0,stop = 6.0,length = 4))
	push!(gs, range(-1,stop =3.0, length=5))
	push!(gs, range(-2,stop = 5.0,length = 7))

	l = Lininterp(vs,gs)
	x = [1.1,3.5,2.8,3.9] 	# == brackets 1,2,4,6
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,4)
	ApproXD.find4DVertices!(l)

	@test length(l.vertex[1]) == 2^l.n
	@test l.vertex[1][1]  == vs[1,2,4,6]  # v0000
	@test l.vertex[1][2]  == vs[2,2,4,6]  # v1000
	@test l.vertex[1][3]  == vs[1,3,4,6]  # v0100
	@test l.vertex[1][4]  == vs[1,2,5,6]  # v0010
	@test l.vertex[1][5]  == vs[1,2,4,7]  # v0001
	@test l.vertex[1][6]  == vs[2,3,4,6]  # v1100
	@test l.vertex[1][7]  == vs[2,2,5,6]  # v1010
	@test l.vertex[1][8]  == vs[2,2,4,7]  # v1001
	@test l.vertex[1][9]  == vs[1,3,5,6]  # v0110
	@test l.vertex[1][10] == vs[1,3,4,7]  # v0101
	@test l.vertex[1][11] == vs[1,2,5,7]  # v0011
	@test l.vertex[1][12] == vs[2,3,5,6]  # v1110
	@test l.vertex[1][13] == vs[2,3,4,7]  # v1101
	@test l.vertex[1][14] == vs[2,2,5,7]  # v1011
	@test l.vertex[1][15] == vs[1,3,5,7]  # v0111
	@test l.vertex[1][16] == vs[2,3,5,7]  # v1111

	x = [3,6,3.0,1.5] 	# == brackets 2,3,4,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,8)
	ApproXD.find4DVertices!(l)

	@test l.vertex[1][1]  == vs[2,3,4,4]  # v0000
	@test l.vertex[1][2]  == vs[3,3,4,4]  # v1000
	@test l.vertex[1][3]  == vs[2,4,4,4]  # v0100
	@test l.vertex[1][4]  == vs[2,3,5,4]  # v0010
	@test l.vertex[1][5]  == vs[2,3,4,5]  # v0001
	@test l.vertex[1][6]  == vs[3,4,4,4]  # v1100
	@test l.vertex[1][7]  == vs[3,3,5,4]  # v1010
	@test l.vertex[1][8]  == vs[3,3,4,5]  # v1001
	@test l.vertex[1][9]  == vs[2,4,5,4]  # v0110
	@test l.vertex[1][10] == vs[2,4,4,5]  # v0101
	@test l.vertex[1][11] == vs[2,3,5,5]  # v0011
	@test l.vertex[1][12] == vs[3,4,5,4]  # v1110
	@test l.vertex[1][13] == vs[3,4,4,5]  # v1101
	@test l.vertex[1][14] == vs[3,3,5,5]  # v1011
	@test l.vertex[1][15] == vs[2,4,5,5]  # v0111
	@test l.vertex[1][16] == vs[3,4,5,5]  # v1111

	# out of bounds
	x = [-1,6.5,-3.0,1.5] 	# == brackets 1,3,1,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (1,11)
	ApproXD.find4DVertices!(l)

	@test l.vertex[1][1]  == vs[1,3,1,4]  # v0000
	@test l.vertex[1][2]  == vs[2,3,1,4]  # v1000
	@test l.vertex[1][3]  == vs[1,4,1,4]  # v0100
	@test l.vertex[1][4]  == vs[1,3,2,4]  # v0010
	@test l.vertex[1][5]  == vs[1,3,1,5]  # v0001
	@test l.vertex[1][6]  == vs[2,4,1,4]  # v1100
	@test l.vertex[1][7]  == vs[2,3,2,4]  # v1010
	@test l.vertex[1][8]  == vs[2,3,1,5]  # v1001
	@test l.vertex[1][9]  == vs[1,4,2,4]  # v0110
	@test l.vertex[1][10] == vs[1,4,1,5]  # v0101
	@test l.vertex[1][11] == vs[1,3,2,5]  # v0011
	@test l.vertex[1][12] == vs[2,4,2,4]  # v1110
	@test l.vertex[1][13] == vs[2,4,1,5]  # v1101
	@test l.vertex[1][14] == vs[2,3,2,5]  # v1011
	@test l.vertex[1][15] == vs[1,4,2,5]  # v0111
	@test l.vertex[1][16] == vs[2,4,2,5]  # v1111
end

@testset "testing find4DVertices!(l) 2 functions" begin
	
	vs1= rand(3,4,5,7)
	vs2= rand(3,4,5,7)
	gs = Array{Float64,1}[]
	push!(gs, range(1.0, stop = 3, length = 3))
	push!(gs, range(2.0,stop = 6.0,length = 4))
	push!(gs, range(-1,stop =3.0, length=5))
	push!(gs, range(-2,stop = 5.0,length = 7))

	l = Lininterp(vs1,vs2,gs)
	x = [1.1,3.5,2.8,3.9] 	# == brackets 1,2,4,6
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,4)
	ApproXD.find4DVertices!(l)

	@test length(l.vertex[1]) == 2^l.n
	@test l.vertex[1][1]  == vs1[1,2,4,6]  # v0000
	@test l.vertex[1][2]  == vs1[2,2,4,6]  # v1000
	@test l.vertex[1][3]  == vs1[1,3,4,6]  # v0100
	@test l.vertex[1][4]  == vs1[1,2,5,6]  # v0010
	@test l.vertex[1][5]  == vs1[1,2,4,7]  # v0001
	@test l.vertex[1][6]  == vs1[2,3,4,6]  # v1100
	@test l.vertex[1][7]  == vs1[2,2,5,6]  # v1010
	@test l.vertex[1][8]  == vs1[2,2,4,7]  # v1001
	@test l.vertex[1][9]  == vs1[1,3,5,6]  # v0110
	@test l.vertex[1][10] == vs1[1,3,4,7]  # v0101
	@test l.vertex[1][11] == vs1[1,2,5,7]  # v0011
	@test l.vertex[1][12] == vs1[2,3,5,6]  # v1110
	@test l.vertex[1][13] == vs1[2,3,4,7]  # v1101
	@test l.vertex[1][14] == vs1[2,2,5,7]  # v1011
	@test l.vertex[1][15] == vs1[1,3,5,7]  # v0111
	@test l.vertex[1][16] == vs1[2,3,5,7]  # v1111

	@test l.vertex[2][1]  == vs2[1,2,4,6]  # v0000
	@test l.vertex[2][2]  == vs2[2,2,4,6]  # v1000
	@test l.vertex[2][3]  == vs2[1,3,4,6]  # v0100
	@test l.vertex[2][4]  == vs2[1,2,5,6]  # v0010
	@test l.vertex[2][5]  == vs2[1,2,4,7]  # v0001
	@test l.vertex[2][6]  == vs2[2,3,4,6]  # v1100
	@test l.vertex[2][7]  == vs2[2,2,5,6]  # v1010
	@test l.vertex[2][8]  == vs2[2,2,4,7]  # v1001
	@test l.vertex[2][9]  == vs2[1,3,5,6]  # v0110
	@test l.vertex[2][10] == vs2[1,3,4,7]  # v0101
	@test l.vertex[2][11] == vs2[1,2,5,7]  # v0011
	@test l.vertex[2][12] == vs2[2,3,5,6]  # v1110
	@test l.vertex[2][13] == vs2[2,3,4,7]  # v1101
	@test l.vertex[2][14] == vs2[2,2,5,7]  # v1011
	@test l.vertex[2][15] == vs2[1,3,5,7]  # v0111
	@test l.vertex[2][16] == vs2[2,3,5,7]  # v1111

	x = [3,6,3.0,1.5] 	# == brackets 2,3,4,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (0,8)
	ApproXD.find4DVertices!(l)

	@test l.vertex[1][1]  == vs1[2,3,4,4]  # v0000
	@test l.vertex[1][2]  == vs1[3,3,4,4]  # v1000
	@test l.vertex[1][3]  == vs1[2,4,4,4]  # v0100
	@test l.vertex[1][4]  == vs1[2,3,5,4]  # v0010
	@test l.vertex[1][5]  == vs1[2,3,4,5]  # v0001
	@test l.vertex[1][6]  == vs1[3,4,4,4]  # v1100
	@test l.vertex[1][7]  == vs1[3,3,5,4]  # v1010
	@test l.vertex[1][8]  == vs1[3,3,4,5]  # v1001
	@test l.vertex[1][9]  == vs1[2,4,5,4]  # v0110
	@test l.vertex[1][10] == vs1[2,4,4,5]  # v0101
	@test l.vertex[1][11] == vs1[2,3,5,5]  # v0011
	@test l.vertex[1][12] == vs1[3,4,5,4]  # v1110
	@test l.vertex[1][13] == vs1[3,4,4,5]  # v1101
	@test l.vertex[1][14] == vs1[3,3,5,5]  # v1011
	@test l.vertex[1][15] == vs1[2,4,5,5]  # v0111
	@test l.vertex[1][16] == vs1[3,4,5,5]  # v1111

	@test l.vertex[2][1]  == vs2[2,3,4,4]  # v0000
	@test l.vertex[2][2]  == vs2[3,3,4,4]  # v1000
	@test l.vertex[2][3]  == vs2[2,4,4,4]  # v0100
	@test l.vertex[2][4]  == vs2[2,3,5,4]  # v0010
	@test l.vertex[2][5]  == vs2[2,3,4,5]  # v0001
	@test l.vertex[2][6]  == vs2[3,4,4,4]  # v1100
	@test l.vertex[2][7]  == vs2[3,3,5,4]  # v1010
	@test l.vertex[2][8]  == vs2[3,3,4,5]  # v1001
	@test l.vertex[2][9]  == vs2[2,4,5,4]  # v0110
	@test l.vertex[2][10] == vs2[2,4,4,5]  # v0101
	@test l.vertex[2][11] == vs2[2,3,5,5]  # v0011
	@test l.vertex[2][12] == vs2[3,4,5,4]  # v1110
	@test l.vertex[2][13] == vs2[3,4,4,5]  # v1101
	@test l.vertex[2][14] == vs2[3,3,5,5]  # v1011
	@test l.vertex[2][15] == vs2[2,4,5,5]  # v0111
	@test l.vertex[2][16] == vs2[3,4,5,5]  # v1111

	# out of bounds
	x = [-1,6.5,-3.0,1.5] 	# == brackets 1,3,1,4
	ApproXD.findBracket!(l,x)
	@test ApproXD.hitmiss(l) == (1,11)
	ApproXD.find4DVertices!(l)

	@test l.vertex[1][1]  == vs1[1,3,1,4]  # v0000
	@test l.vertex[1][2]  == vs1[2,3,1,4]  # v1000
	@test l.vertex[1][3]  == vs1[1,4,1,4]  # v0100
	@test l.vertex[1][4]  == vs1[1,3,2,4]  # v0010
	@test l.vertex[1][5]  == vs1[1,3,1,5]  # v0001
	@test l.vertex[1][6]  == vs1[2,4,1,4]  # v1100
	@test l.vertex[1][7]  == vs1[2,3,2,4]  # v1010
	@test l.vertex[1][8]  == vs1[2,3,1,5]  # v1001
	@test l.vertex[1][9]  == vs1[1,4,2,4]  # v0110
	@test l.vertex[1][10] == vs1[1,4,1,5]  # v0101
	@test l.vertex[1][11] == vs1[1,3,2,5]  # v0011
	@test l.vertex[1][12] == vs1[2,4,2,4]  # v1110
	@test l.vertex[1][13] == vs1[2,4,1,5]  # v1101
	@test l.vertex[1][14] == vs1[2,3,2,5]  # v1011
	@test l.vertex[1][15] == vs1[1,4,2,5]  # v0111
	@test l.vertex[1][16] == vs1[2,4,2,5]  # v1111

	@test l.vertex[2][1]  == vs2[1,3,1,4]  # v0000
	@test l.vertex[2][2]  == vs2[2,3,1,4]  # v1000
	@test l.vertex[2][3]  == vs2[1,4,1,4]  # v0100
	@test l.vertex[2][4]  == vs2[1,3,2,4]  # v0010
	@test l.vertex[2][5]  == vs2[1,3,1,5]  # v0001
	@test l.vertex[2][6]  == vs2[2,4,1,4]  # v1100
	@test l.vertex[2][7]  == vs2[2,3,2,4]  # v1010
	@test l.vertex[2][8]  == vs2[2,3,1,5]  # v1001
	@test l.vertex[2][9]  == vs2[1,4,2,4]  # v0110
	@test l.vertex[2][10] == vs2[1,4,1,5]  # v0101
	@test l.vertex[2][11] == vs2[1,3,2,5]  # v0011
	@test l.vertex[2][12] == vs2[2,4,2,4]  # v1110
	@test l.vertex[2][13] == vs2[2,4,1,5]  # v1101
	@test l.vertex[2][14] == vs2[2,3,2,5]  # v1011
	@test l.vertex[2][15] == vs2[1,4,2,5]  # v0111
	@test l.vertex[2][16] == vs2[2,4,2,5]  # v1111
end

@testset "testing getValue 1D" begin

	lbs = [1.0]
	ubs = [3.0]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 4))

	myfun(i1) = 0.75*i1 
	
	vs = Float64[ myfun(i) for i in gs[1]]

	l = Lininterp(vs,gs)

	# check value on bounds
	@test getValue(l,lbs)[1] == vs[1]
	@test getValue(l,ubs)[1] == vs[4]
	@test getValue(l,[5.0])[1] == vs[4]
	@test getValue(l,[1.0])[1]== vs[1]

	# check values out of bounds
	@test getValue(l,[-1.0])[1] == vs[1]
	@test getValue(l,[200.0])[1] == vs[4]
	println(l)

	# close to bounds
	y=2.9
	@test isapprox(getValue(l,[y])[1] - myfun(y) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[y])[1] - myfun(y) ,0.0,atol=1e-6)


	# check at random vals in interval
	lbs = [1.0]
	ubs = [33.0]
	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1], length = 40))

	vs = Float64[ myfun(i) for i in gs[1]]

	l = Lininterp(vs,gs)
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		@test isapprox(getValue(l,[x])[1] ,myfun(x),atol=1e-6)
	end
	println(l)
end

@testset "testing getValue! 1D" begin

	myfun(i1) = 0.75*i1 
	
	lbs = [1.0]
	ubs = [33.0]
	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1], length = 40))

	vs = Float64[ myfun(i) for i in gs[1]]

	l = Lininterp(vs,gs)
	y = [0.0]
	# check at increasing vals in interval
	x = range(lbs[1],stop = ubs[1],length = 300)
	for i in 1:length(x)
		getValue!(y,l,[x[i]],[1])
		@test isapprox(y[1] ,myfun(x[i]),atol=1e-8)
	end
	println(l)
	resetCache!(l)
	# check at random vals in interval
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		getValue!(y,l,[x],[1])
		@test isapprox(y[1] ,myfun(x),atol=1e-8)
	end
	println(l)
end



@testset "testing getValue! 1D on irregularly spaced grid" begin

	lb = 1.0
	ub = 3.0
	n = 40

	xg = zeros(n)
	xg[1] = log(lb + 1) 
	xg[n] = log(ub + 1) 
	xg    = range(xg[1],stop = xg[n],length = n)
	xg    = exp.(xg) .- 1  

	gs = Array{Float64,1}[]
	push!(gs, xg)

	# I'm using a linear function here: the more non-linear your function,
	# the worse a linear approximation will perform - this is obvious.
	myfun(i1) = 0.75*i1 
	

	vs = Float64[ myfun(i) for i in gs[1]]

	l = Lininterp(vs,gs)
	y = [0.0]
	# check at increasing vals in interval
	x = range(lb,stop = ub,length = 300)
	for i in 1:length(x)
		getValue!(y,l,[x[i]],[1])
		@test isapprox(y[1] ,myfun(x[i]),atol=1e-8)
		# println(y[1] - myfun(x[i]))
	end
	println(l)
	resetCache!(l)
	# check at random vals in interval
	for i in 1:30
		x = rand() * (ub-lb) + lb
		getValue!(y,l,[x],[1])
		@test isapprox(y[1] ,myfun(x),atol=1e-8)
	end
	println(l)
end

@testset "testing getValue2D" begin

	lbs = [1.0,2.0]
	ubs = [3.0,5.0]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 3))
	push!(gs, range(lbs[2],stop = ubs[2],length = 4))

	myfun(i1,i2) = 0.5*i1 + 2*i2 
	
	vs = Float64[ myfun(i,j) for i in gs[1], j in gs[2]]

	l = Lininterp(vs,gs)

	# check value on bounds
	@test getValue(l,lbs)[1] == vs[1,1]
	@test getValue(l,ubs)[1] == vs[3,4]
	@test getValue(l,[1.0,5.0])[1] == vs[1,4]
	@test getValue(l,[1.0,5.0])[1]== vs[1,4]

	# check values out of bounds
	@test getValue(l,[-1.0,2.0])[1] == vs[1,1]
	@test getValue(l,[1.0,200.0])[1] == vs[1,4]
	println(l)

	# close to bounds
	x=1.99
	y=4.9
	@test isapprox(getValue(l,[x,y])[1] - myfun(x,y) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y])[1] - myfun(x,y) ,0.0,atol=1e-6)

	x=1.4861407584066377
	y=3.5646251730324234
	@test isapprox(getValue(l,[x,y])[1] - myfun(x,y) ,0.0,atol=1e-6)

	x = 2.53
	y = 2.58
	@test isapprox(getValue(l,[x,y])[1] - myfun(x,y) ,0.0,atol=1e-6)



	# check at random vals in interval
	lbs = [1.0,2.0]
	ubs = [33.0,5.0]
	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 30))
	push!(gs, range(lbs[2],stop = ubs[2],length = 40))

	vs = Float64[ myfun(i,j) for i in gs[1], j in gs[2]]

	l = Lininterp(vs,gs)
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		y = rand() * (ubs[2]-lbs[2]) + lbs[2]
		@test isapprox(getValue(l,[x,y])[1] ,myfun(x,y),atol=1e-6)
	end
	println(l)
end




@testset "testing getValue3D" begin

	lbs = [1.0,2.0,-1]
	ubs = [3.0,5.0,3]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 3))
	push!(gs, range(lbs[2],stop = ubs[2],length = 4))
	push!(gs, range(lbs[3],stop = ubs[3],length = 5))

	myfun(i1,i2,i3) = 0.5*i1 + 2*i2 + 3*i3
	
	vs = Float64[ myfun(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs,gs)

	# check value on bounds
	@test getValue(l,lbs)[1] == vs[1,1,1]
	@test getValue(l,ubs)[1] == vs[3,4,5]
	@test getValue(l,[1.0,5.0,-1])[1] == vs[1,4,1]
	@test getValue(l,[1.0,5.0,3])[1] == vs[1,4,5]

	# check values out of bounds
	@test getValue(l,[-1.0,2.0,-1])[1] == vs[1,1,1]
	@test getValue(l,[1.0,200.0,-1])[1] == vs[1,4,1]
	println(l)
	# ApproXD.resetCache!(l)

	# close to bounds
	x=1.99
	y=4.9
	z=2.9
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)

	x=1.4861407584066377
	y=3.5646251730324234
	z=2.7832610606532944
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)

	x = 2.53
	y = 2.58
	z = -0.87
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)

	# check at random vals in interval
	lbs = [1.0,2.0,-10]
	ubs = [33.0,5.0,3]
	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 30))
	push!(gs, range(lbs[2],stop = ubs[2],length = 40))
	push!(gs, range(lbs[3],stop = ubs[3],length = 50))

	vs = Float64[ myfun(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs,gs)
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		y = rand() * (ubs[2]-lbs[2]) + lbs[2]
		z = rand() * (ubs[3]-lbs[3]) + lbs[3]
		@test isapprox(getValue(l,[x,y,z])[1] ,myfun(x,y,z),atol=1e-6)
	end
	println(l)
end


@testset "testing getValue3D on 2 functions" begin

	lbs = [1.0,2.0,-1]
	ubs = [3.0,5.0,3]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 3))
	push!(gs, range(lbs[2],stop = ubs[2],length = 4))
	push!(gs, range(lbs[3],stop = ubs[3],length = 5))

	myfun1(i1,i2,i3) = 0.5*i1 + 2*i2 + 3*i3
	myfun2(i1,i2,i3) = 3*i1 + pi*i2 + 0.1*i3
	
	vs1 = Float64[ myfun1(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]
	vs2 = Float64[ myfun2(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs1,vs2,gs)

	# check value on bounds
	@test getValue(l,lbs) .- [vs1[1,1,1],vs2[1,1,1]] == zeros(2)
	@test getValue(l,ubs) .- [vs1[3,4,5],vs2[3,4,5]] == zeros(2)
	@test getValue(l,[1.0,5.0,-1])[1] == vs1[1,4,1]
	@test getValue(l,[1.0,5.0,-1])[2] == vs2[1,4,1]
	@test getValue(l,[1.0,5.0,3])[1] == vs1[1,4,5]
	@test getValue(l,[1.0,5.0,3])[2] == vs2[1,4,5]

	# check values out of bounds
	@test getValue(l,[-1.0,2.0,-1])[1] == vs1[1,1,1]
	@test getValue(l,[-1.0,2.0,-1])[2] == vs2[1,1,1]
	@test getValue(l,[1.0,200.0,-1])[1] == vs1[1,4,1]
	@test getValue(l,[1.0,200.0,-1])[2] == vs2[1,4,1]
	println(l)
	# ApproXD.resetCache!(l)

	# close to bounds
	x=1.99
	y=4.9
	z=2.9
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)

	x=1.4861407584066377
	y=3.5646251730324234
	z=2.7832610606532944
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)

	x = 2.53
	y = 2.58
	z = -0.87
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)

	# check at random vals in interval
	lbs = [1.0,2.0,-10]
	ubs = [33.0,5.0,3]
	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 30))
	push!(gs, range(lbs[2],stop = ubs[2],length = 40))
	push!(gs, range(lbs[3],stop = ubs[3],length = 50))

	vs1 = Float64[ myfun1(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]
	vs2 = Float64[ myfun2(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs1,vs2,gs)
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		y = rand() * (ubs[2]-lbs[2]) + lbs[2]
		z = rand() * (ubs[3]-lbs[3]) + lbs[3]
		@test isapprox(getValue(l,[x,y,z])[1] ,myfun1(x,y,z),atol=1e-6)
		@test isapprox(getValue(l,[x,y,z])[2] ,myfun2(x,y,z),atol=1e-6)
	end
	println(l)
end

@testset "testing getValue! 3D on 2 functions with ifunc switch" begin

	lbs = [1.0,2.0,-1]
	ubs = [3.0,5.0,3]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1],length = 3))
	push!(gs, range(lbs[2],stop = ubs[2],length = 4))
	push!(gs, range(lbs[3],stop = ubs[3],length = 5))

	myfun1(i1,i2,i3) = 0.5*i1 + 2*i2 + 3*i3
	myfun2(i1,i2,i3) = 3*i1 + pi*i2 + 0.1*i3
	
	vs1 = Float64[ myfun1(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]
	vs2 = Float64[ myfun2(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs1,vs2,gs)

	@test_throws MethodError getValue!(l,lbs,[1,2,3])

	# check value on bounds
	y = [0.0]
	getValue!(y,l,lbs,[1])
	@test y == [vs1[1,1,1]]

	getValue!(y,l,lbs,[2])
	@test y == [vs2[1,1,1]]
	getValue!(y,l,ubs,[1])
	@test y == [vs1[3,4,5]]
	getValue!(y,l,ubs,[2]) 
	@test y == [vs2[3,4,5]]

	getValue!(y,l,[1.0,5.0,-1],[1])
	@test y == [vs1[1,4,1]]
	getValue!(y,l,[1.0,5.0,-1],[2])
	@test y == [vs2[1,4,1]]
	getValue!(y,l,[1.0,5.0,3],[1]) 
	@test y  == [vs1[1,4,5]]
	getValue!(y,l,[1.0,5.0,3],[2])
	@test y  == [vs2[1,4,5]]

	# check values out of bounds
	getValue!(y,l,[-1.0,2.0,-1],[1])
	@test y == [vs1[1,1,1]]
	getValue!(y,l,[-1.0,2.0,-1],[2]) 
	@test y == [vs2[1,1,1]]
	y2 = [0.0,0.0]
	getValue!(y2,l,[1.0,200.0,-1],[1,2])
	@test y2 == [vs1[1,4,1],vs2[1,4,1]]
	getValue!(y,l,[1.0,200.0,-1],[2])
	@test y == [vs2[1,4,1]]
	println(l)
	# ApproXD.resetCache!(l)

	# close to bounds
	x=1.99
	y=4.9
	z=2.9
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)

	x=1.4861407584066377
	y=3.5646251730324234
	z=2.7832610606532944
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)

	x = 2.53
	y = 2.58
	z = -0.87
	@test isapprox(getValue(l,[x,y,z])[1] - myfun1(x,y,z) ,0.0,atol=1e-6)
	@test isapprox(getValue(l,[x,y,z])[2] - myfun2(x,y,z) ,0.0,atol=1e-6)

	# check at random vals in interval
	lbs = [1.0,2.0,-10]
	ubs = [33.0,5.0,3]
	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1], length = 30))
	push!(gs, range(lbs[2],stop = ubs[2], length = 40))
	push!(gs, range(lbs[3],stop = ubs[3], length = 50))

	vs1 = Float64[ myfun1(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]
	vs2 = Float64[ myfun2(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs1,vs2,gs)
	yout = [0.0,0.0]
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		y = rand() * (ubs[2]-lbs[2]) + lbs[2]
		z = rand() * (ubs[3]-lbs[3]) + lbs[3]
		getValue!(yout,l,[x,y,z],[1,2])
		@test isapprox(yout[1] ,myfun1(x,y,z),atol=1e-6)
		@test isapprox(yout[2] ,myfun2(x,y,z),atol=1e-6)
	end
	println(l)
end

@testset "testing getValue hit/miss for 3D" begin

	lbs = [1.0,2.0,-1]
	ubs = [3.0,5.0,3]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1], length = 3))
	push!(gs, range(lbs[2],stop = ubs[2], length = 4))
	push!(gs, range(lbs[3],stop = ubs[3], length = 5))

	myfun(i1,i2,i3) = i1 + 2*i2 + 3*i3
	
	vs = Float64[ myfun(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = Lininterp(vs,gs)
	@test hitmiss(l) == (0,0)
	
	x=1.99
	y=4.9
	z=2.9
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)
	@test hitmiss(l) == (0,3)
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)
	@test hitmiss(l) == (3,3)
	x=2.99
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)
	@test hitmiss(l) == (5,4)
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)
	@test hitmiss(l) == (8,4)
	y=2.11
	@test isapprox(getValue(l,[x,y,z])[1] - myfun(x,y,z) ,0.0,atol=1e-6)
	@test hitmiss(l) == (10,5)

end



@testset "testing getValue4D on 3 functions" begin

	lbs = [1.0,2.0,-1,4]
	ubs = [3.0,5.0,3,18.0]

	gs = Array{Float64,1}[]
	push!(gs, range(lbs[1],stop = ubs[1], length = 3))
	push!(gs, range(lbs[2],stop = ubs[2], length = 4))
	push!(gs, range(lbs[3],stop = ubs[3], length = 5))
	push!(gs, range(lbs[4],stop = ubs[4], length = 9))

	myfun1(i1,i2,i3,i4) = 0.5*i1 + 2*i2 + 3*i3 + i4/0.98
	myfun2(i1,i2,i3,i4) = i1 + 0.2*i2 + 3*i3 + i4
	myfun3(i1,i2,i3,i4) = 0.7*i1 + 1.5*i2 + 2*i3 + i4/0.2
	
	vs1 = Float64[ myfun1(i,j,k,m) for i in gs[1], j in gs[2], k in gs[3], m in gs[4] ]
	vs2 = Float64[ myfun2(i,j,k,m) for i in gs[1], j in gs[2], k in gs[3], m in gs[4] ]
	vs3 = Float64[ myfun3(i,j,k,m) for i in gs[1], j in gs[2], k in gs[3], m in gs[4] ]

		vs = Array{Float64}[]
		push!(vs,vs1)
		push!(vs,vs2)
		push!(vs,vs3)

		l = Lininterp(vs,gs)

	@testset "returning all functions stored on L" begin

		# check value on bounds
		v = getValue(l,lbs)
		@test v[1] == vs1[1,1,1,1]
		@test v[2] == vs2[1,1,1,1]
		@test v[3] == vs3[1,1,1,1]

		v = getValue(l,ubs)
		@test v[1] == vs1[3,4,5,9]
		@test v[2] == vs2[3,4,5,9]
		@test v[3] == vs3[3,4,5,9]

		v = getValue(l,[1.0,5.0,-1,4])
		@test v[1] == vs1[1,4,1,1]
		@test v[2] == vs2[1,4,1,1]
		@test v[3] == vs3[1,4,1,1]

		v = getValue(l,[1.0,5.0,3,18])
		@test v[1] == vs1[1,4,5,9]
		@test v[2] == vs2[1,4,5,9]
		@test v[3] == vs3[1,4,5,9]

		# check values out of bounds
		v = getValue(l,[-1.0,2.0,-1,4])
		@test v[1] == vs1[1,1,1,1]
		@test v[2] == vs2[1,1,1,1]
		@test v[3] == vs3[1,1,1,1]

		v = getValue(l,[1.0,200.0,-1,18])
		@test v[1] == vs1[1,4,1,9]
		@test v[2] == vs2[1,4,1,9]
		@test v[3] == vs3[1,4,1,9]
		println(l)
		# ApproXD.resetCache!(l)

		# close to bounds
		x=1.99
		y=4.9
		z=2.9
		w=17.95
		v = getValue(l,[x,y,z,w])
		@test isapprox(v[1] - myfun1(x,y,z,w) ,0.0,atol=1e-6)
		@test isapprox(v[2] - myfun2(x,y,z,w) ,0.0,atol=1e-6)
		@test isapprox(v[3] - myfun3(x,y,z,w) ,0.0,atol=1e-6)

		x=1.4861407584066377
		y=3.5646251730324234
		z=2.7832610606532944
		w=10.05
		v = getValue(l,[x,y,z,w])
		@test isapprox(v[1] - myfun1(x,y,z,w) ,0.0,atol=1e-6)
		@test isapprox(v[2] - myfun2(x,y,z,w) ,0.0,atol=1e-6)
		@test isapprox(v[3] - myfun3(x,y,z,w) ,0.0,atol=1e-6)

		x = 2.53
		y = 2.58
		z = -0.87
		w = 5.77
		v = getValue(l,[x,y,z,w])
		@test isapprox(v[1] - myfun1(x,y,z,w) ,0.0,atol=1e-6)
		@test isapprox(v[2] - myfun2(x,y,z,w) ,0.0,atol=1e-6)
		@test isapprox(v[3] - myfun3(x,y,z,w) ,0.0,atol=1e-6)

		# check at random vals in interval

		for i in 1:100
			x = rand() * (ubs[1]-lbs[1]) + lbs[1]
			y = rand() * (ubs[2]-lbs[2]) + lbs[2]
			z = rand() * (ubs[3]-lbs[3]) + lbs[3]
			w = rand() * (ubs[4]-lbs[4]) + lbs[4]
			v = getValue(l,[x,y,z,w])
			@test isapprox(v[1] - myfun1(x,y,z,w) ,0.0,atol=1e-6)
			@test isapprox(v[2] - myfun2(x,y,z,w) ,0.0,atol=1e-6)
			@test isapprox(v[3] - myfun3(x,y,z,w) ,0.0,atol=1e-6)
		end
	end

	@testset "returning selected functions stored on L" begin
		f = [myfun1;myfun2;myfun3]
		for i in 1:10
			x = rand() * (ubs[1]-lbs[1]) + lbs[1]
			y = rand() * (ubs[2]-lbs[2]) + lbs[2]
			z = rand() * (ubs[3]-lbs[3]) + lbs[3]
			w = rand() * (ubs[4]-lbs[4]) + lbs[4]
			gets = ApproXD.sample(1:l.nfunc,rand(1:l.nfunc))
			v = zeros(length(gets))
			getValue!(v,l,[x,y,z,w],gets)
			for ig in 1:length(gets)
				@test isapprox(v[ig] - f[gets[ig]](x,y,z,w) ,0.0,atol=1e-6)
			end
		end
	end

end


end