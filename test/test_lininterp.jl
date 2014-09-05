

module test_lininterp

using ApproXD, FactCheck

facts("constructor for lininterp") do

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, linspace(1.0,3,3))
	push!(gs, linspace(2.0,3,4))
	push!(gs, linspace(-1,3.0,5))

	l = lininterp(vs,gs)

	@fact ApproXD.getDims(l) => [3,4,5]
	@fact l.n => 3
	@fact l.infs => zeros(3)
	@fact l.sups => zeros(3)
	@fact l.hits => 0
	@fact l.miss => 0
	@fact l.hitnow => false
	@fact l.vals => vs
	@fact l.grids => gs
	@fact ApproXD.getCache(l) => [1,1,1]

	# errors
	pop!(gs)   # missing grid
	@fact_throws lininterp(vs,gs)   
	push!(gs, linspace(-1,3.0,4)) 	# wrong length
	@fact_throws lininterp(vs,gs)
	pop!(gs)
	push!(gs, [1.1,0,2,3,4]) 	# not sorted
	@fact_throws lininterp(vs,gs)

end

facts("testing Cache functions") do

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, linspace(1.0,3,3))
	push!(gs, linspace(2.0,3,4))
	push!(gs, linspace(-1,3.0,5))

	l = lininterp(vs,gs)
	myc = [2,1,3]
	l.cache = myc	# manually set a cache

	for i in 1:length(gs)
		@fact ApproXD.getCache(l,i) => myc[i]
		@fact ApproXD.getCachedVal(l,i) => gs[i][myc[i]]
		@fact ApproXD.getNextCachedVal(l,i) => gs[i][myc[i]+1]
	end


end

facts("testing findBracket!(x)") do

	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, linspace(1.0,3,3))
	push!(gs, linspace(2.0,6.0,4))
	push!(gs, linspace(-1,3.0,5))

	l = lininterp(vs,gs)
	@fact l.hascache => false

	x = [1.1,3.5,2.8] 	# => brackets 1,2,4
	ApproXD.findBracket!(l,x)

	# brackets is index of lower bound of each bracket
	brackets = [1,2,4]

	@fact ApproXD.getCache(l) => brackets
	@fact ApproXD.hitmiss(l) => (0,3)
	for i in 1:3
		@fact l.infs[i] => gs[i][brackets[i]]
		@fact l.sups[i] => gs[i][brackets[i]+1]
	end

	@fact l.hascache => true
	ApproXD.findBracket!(l,x)
	@fact ApproXD.getCache(l) => brackets
	@fact ApproXD.hitmiss(l) => (3,3)

	# value x1 and x2 out of bounds
	x = [0.9,7.5,2.8] 	# => brackets 1,3,4
	ApproXD.findBracket!(l,x)

	brackets = [1,3,4]
	@fact ApproXD.getCache(l) => brackets
	@fact ApproXD.hitmiss(l) => (4,5)
	@fact l.infs[1] => gs[1][1]
	@fact l.infs[2] => gs[2][end-1]
	@fact l.infs[3] => gs[3][brackets[3]]
	@fact l.sups[1] => gs[1][2]
	@fact l.sups[2] => gs[2][end]
	@fact l.sups[3] => gs[3][brackets[3]+1]

	# side-effect: x is forced into bounds!
	# may or may not be desirable.
	@fact x[1] => l.infs[1]
	@fact x[2] => l.sups[2]
end

facts("testing find3DVertices!(l)") do
	
	vs = rand(3,4,5)
	gs = Array{Float64,1}[]
	push!(gs, linspace(1.0,3,3))
	push!(gs, linspace(2.0,6.0,4))
	push!(gs, linspace(-1,3.0,5))

	l = lininterp(vs,gs)
	x = [1.1,3.5,2.8] 	# => brackets 1,2,4
	ApproXD.findBracket!(l,x)
	@fact ApproXD.hitmiss(l) => (0,3)
	ApproXD.find3DVertices!(l)

	@fact length(l.vertex) => 8
	@fact l.vertex[1] => vs[1,2,4]  # v000   
	@fact l.vertex[2] => vs[2,2,4]  # v100
	@fact l.vertex[3] => vs[1,3,4]  # v010
	@fact l.vertex[4] => vs[1,2,5]  # v001
	@fact l.vertex[5] => vs[2,3,4]  # v110
	@fact l.vertex[6] => vs[2,2,5]  # v101
	@fact l.vertex[7] => vs[1,3,5]  # v011
	@fact l.vertex[8] => vs[2,3,5]  # v111

	x = [3,6,3.0] 	# => brackets 2,3,4
	ApproXD.findBracket!(l,x)
	@fact ApproXD.hitmiss(l) => (0,6)
	@fact l.hitnow => false
	ApproXD.find3DVertices!(l)

	@fact l.vertex[1] => vs[2,3,4]  # v000   
	@fact l.vertex[2] => vs[3,3,4]  # v100
	@fact l.vertex[3] => vs[2,4,4]  # v010
	@fact l.vertex[4] => vs[2,3,5]  # v001
	@fact l.vertex[5] => vs[3,4,4]  # v110
	@fact l.vertex[6] => vs[3,3,5]  # v101
	@fact l.vertex[7] => vs[2,4,5]  # v011
	@fact l.vertex[8] => vs[3,4,5]  # v111

	# out of bounds
	x = [-1,6.5,-3.0] 	# => brackets 1,3,1
	ApproXD.findBracket!(l,x)
	@fact ApproXD.hitmiss(l) => (0,9)
	@fact l.hitnow => false
	ApproXD.find3DVertices!(l)

	@fact l.vertex[1] => vs[1,3,1]  # v000   
	@fact l.vertex[2] => vs[2,3,1]  # v100
	@fact l.vertex[3] => vs[1,4,1]  # v010
	@fact l.vertex[4] => vs[1,3,2]  # v001
	@fact l.vertex[5] => vs[2,4,1]  # v110
	@fact l.vertex[6] => vs[2,3,2]  # v101
	@fact l.vertex[7] => vs[1,4,2]  # v011
	@fact l.vertex[8] => vs[2,4,2]  # v111
end


facts("testing eval3D") do

	lbs = [1.0,2.0,-1]
	ubs = [3.0,5.0,3]

	gs = Array{Float64,1}[]
	push!(gs, linspace(lbs[1],ubs[1],3))
	push!(gs, linspace(lbs[2],ubs[2],4))
	push!(gs, linspace(lbs[3],ubs[3],5))

	myfun(i1,i2,i3) = i1 + 2*i2 + 3*i3
	
	vs = Float64[ myfun(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = lininterp(vs,gs)

	# check value on bounds
	@fact eval3D(l,lbs) => vs[1,1,1]
	@fact eval3D(l,ubs) => vs[3,4,5]
	@fact eval3D(l,[1.0,5.0,-1]) => vs[1,4,1]
	@fact eval3D(l,[1.0,5.0,3]) => vs[1,4,5]

	# check values out of bounds
	@fact eval3D(l,[-1.0,2.0,-1]) => vs[1,1,1]
	println("hitnow = $(l.hitnow)")
	println("cache = $(l.cache)")
	@fact eval3D(l,[1.0,200.0,-1]) => vs[1,4,1]
	println(l)
	println("cache = $(l.cache)")
	println(l.vertex)
	println("hitnow = $(l.hitnow)")
	println("hascache = $(l.hascache)")
	# ApproXD.resetCache!(l)

	# close to bounds
	x=1.99
	y=4.9
	z=2.9
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)

	x=1.4861407584066377
	y=3.5646251730324234
	z=2.7832610606532944
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)

	x = 2.53
	y = 2.58
	z = -0.87
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)

	# check at random vals in interval
	lbs = [1.0,2.0,-10]
	ubs = [33.0,5.0,3]
	gs = Array{Float64,1}[]
	push!(gs, linspace(lbs[1],ubs[1],30))
	push!(gs, linspace(lbs[2],ubs[2],40))
	push!(gs, linspace(lbs[3],ubs[3],50))

	vs = Float64[ myfun(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = lininterp(vs,gs)
	for i in 1:30
		x = rand() * (ubs[1]-lbs[1]) + lbs[1]
		y = rand() * (ubs[2]-lbs[2]) + lbs[2]
		z = rand() * (ubs[3]-lbs[3]) + lbs[3]
		@fact eval3D(l,[x,y,z]) => roughly(myfun(x,y,z),atol=1e-6)
	end
	println(l)
end

facts("testing eval3D hit/miss") do

	lbs = [1.0,2.0,-1]
	ubs = [3.0,5.0,3]

	gs = Array{Float64,1}[]
	push!(gs, linspace(lbs[1],ubs[1],3))
	push!(gs, linspace(lbs[2],ubs[2],4))
	push!(gs, linspace(lbs[3],ubs[3],5))

	myfun(i1,i2,i3) = i1 + 2*i2 + 3*i3
	
	vs = Float64[ myfun(i,j,k) for i in gs[1], j in gs[2], k in gs[3] ]

	l = lininterp(vs,gs)
	@fact hitmiss(l) => (0,0)
	
	x=1.99
	y=4.9
	z=2.9
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)
	@fact hitmiss(l) => (0,3)
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)
	@fact hitmiss(l) => (3,3)
	x=2.99
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)
	@fact hitmiss(l) => (5,4)
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)
	@fact hitmiss(l) => (8,4)
	y=2.11
	@fact eval3D(l,[x,y,z]) - myfun(x,y,z) => roughly(0.0,atol=1e-6)
	@fact hitmiss(l) => (10,5)

end



end
