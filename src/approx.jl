





# approximating tensor grids
# computing the approximating coefficients on 
# a tensor product of basis functions
# NOTE: the fastest varying index in v is the one with highest 
# index in ibm
function getTensorCoef(ibm::Dict{Integer, Array{T, 2}}, v::Vector{T}) where T <: Real

	# ibm are usually inverse basis matrices

	nall = length(v)	# length value vector
	nbm  = length(ibm)	# number of basis functions / dimensions of v

	# dim(matrices) must be same as length v
	prod = 1
	for (k,va) in ibm
		prod *= size(va,1)
		if !(size(va,1) == size(va,2))
			throw(ArgumentError("all matrices must be square"))
		end
	end
	if !(prod == nall)
		throw(ArgumentError("sizes of matrices not consistent with v"))
	end

	vtmp = 0.0

	# compute product for first matrix
	ks = sort(collect(keys(ibm)))
	v0    = copy(v)
	v1    = zeros(nall)
	stemp = ibm[ks[1]]
	n     = size(stemp,1)
	m     = round(Int,nall / n)
	for i in 1:m, j=1:n
		vtmp = 0.0
		for ji=1:n
			@inbounds vtmp += stemp[j,ji] * v0[n*(i-1) + ji]
		end
		v1[m*(j-1) + i] = vtmp
		# v[m*(j-1) + i] = vtmp
	end

	# compute for all other matrics
	if nbm > 1
		for imat in 2:nbm
			v0    = copy(v1)
			# v0    = copy(v)
			stemp = ibm[ks[imat]]
			n     = size(stemp,1)
			m     = round(Int,nall / n)
			# fill!(v1,0.0)
			# for i in 1:m, j=1:n
			for i in 1:m, j=1:n
				vtmp = 0.0
				for ji=1:n
					@inbounds vtmp += stemp[j,ji] * v0[n*(i-1) + ji]
				end
				# v[m*(j-1) + i] = vtmp
				v1[m*(j-1) + i] = vtmp
			end
		end
	end
	return v1
end



function evalTensor2(mat1::Vector{T}, mat2::Vector{T}, c::Vector) where T

	# TODO
	# sparse: if you had a type BSpline with field "nonzero" that
	# gives the index of the nonzero basis functions, this could
	# be sped up dramatically

	n1 = size(mat1,1)
	n2 = size(mat2,1)
	m1 = size(mat1,2)
	m2 = size(mat2,2)
	row_offset = 0
	col_offset = 0
	factor = 0.0

	r = zeros(n1 * n2)

	# loop over rows of 1
	for row_idx1 in 1:n1
		row_offset = (row_idx1-1) * n2

		# loop over rows of 2
		for row_idx2 in 1:n2

			# loop over cols of 1
			for col_idx1 in 1:m1
				col_offset = (col_idx1-1) * m2
				factor = mat1[row_idx1,col_idx1]

				# loop over cols of 2
				for col_idx2 in 1:m2

					# println("row_offset = $row_offset")
					# println("row_idx1 = $row_idx1")
					# println("row_idx2 = $row_idx2")
					# println("row_offset + row_idx1 = $(row_offset + row_idx1)")
					# println("factor = $factor")

					r[row_offset + row_idx2] += factor * mat2[row_idx2,col_idx2] * c[col_offset + col_idx2]
				end
			end
		end
	end
	return r
end


function evalTensor2(mat1::SparseMatrixCSC{T, Int64}, mat2::SparseMatrixCSC{T, Int64}, c::Vector{T}) where T

	# assume that the basis matrix is stored colwise
	# i.e. the first column are the basis functions
	# evaluated at the first point.

	# FIXME
	# julia 0.4 will introduce sparse row compression
	# so this can be changed in BSplines.jl

	if size(mat1)[2] > 1
		error("only doing column vector so far, sorry")
	end

	m1 = length(mat1)
	m2 = length(mat2)
	offset = 0
	factor = 0.0

	r = 0.0

	# full idiom is 
	# for c in 1:size(mat1)[2]
	#    non_zero_indices[c] = mat1[mat1.colptr[c] : mat1.colptr[c+1]-1]
	#

	# for now (one column only) this is safe:
	# non_zero_indices = mat.rowval

	# loop over non-zero entries of mat1
	for idx1 in mat1.rowval
		offset = (idx1-1) * m2
		factor = mat1[idx1]

		# loop over non-zeros of mat2
		for idx2 in mat2.rowval

			r += factor * mat2[idx2] * c[offset + idx2]
		end
	end
	return r
end


function evalTensor3(mat1::SparseMatrixCSC{T, Int64}, mat2::SparseMatrixCSC{T, Int64}, mat3::SparseMatrixCSC{T, Int64}, c::Vector{T}) where T

	n1 = size(mat1,1)
	n2 = size(mat2,1)
	n3 = size(mat3,1)
	m1 = size(mat1,2)
	m2 = size(mat2,2)
	m3 = size(mat3,2)
	row_offset1= 0
	col_offset1= 0
	factor1= 0.0
	row_offset2= 0
	col_offset2= 0
	factor2= 0.0

	r = zeros(n1 * n2 * n3)

	# loop over rows of 1
	for row_idx1 in 1:n1
		row_offset1 = (row_idx1-1) * n2

		# loop over rows of 2
		for row_idx2 in 1:n2
			row_offset2 = (row_offset1 + (row_idx2-1)) * n3

			# loop rows 3
			for row_idx3 in 1:n3
				# println(row_offset2 + row_idx3)

				# loop over cols of 1
				for col_idx1 in 1:m1
					col_offset1 = (col_idx1-1) * m2
					factor1 = mat1[row_idx1,col_idx1]

					# loop over cols of 2
					for col_idx2 in 1:m2
						col_offset2 = (col_offset1 + (col_idx2-1)) * m3
						factor2 = factor1 * mat2[row_idx2,col_idx2]

						# cols mat3
						for col_idx3 in 1:m3

							r[row_offset2 + row_idx3] += factor2 * mat3[row_idx3,col_idx3] * c[col_offset2 + col_idx3]
						end
					end
				end
			end
		end
	end
	return r
end


function evalTensor3(mat1::Array{Float64,1},mat2::Array{Float64,1},mat3::Array{Float64,1},c::Vector)

	m1 = length(mat1)
	m2 = length(mat2)
	m3 = length(mat3)
	col_offset1= 0
	factor1= 0.0
	col_offset2= 0
	factor2= 0.0

	r = 0.0

	# loop over cols of 1
	for col_idx1 in 1:m1
		col_offset1 = (col_idx1-1) * m2
		factor1 = mat1[col_idx1]

		# loop over cols of 2
		for col_idx2 in 1:m2
			col_offset2 = (col_offset1 + (col_idx2-1)) * m3
			factor2 = factor1 * mat2[col_idx2]

			# cols mat3
			for col_idx3 in 1:m3
				r += factor2 * mat3[col_idx3] * c[col_offset2 + col_idx3]
			end
		end
	end
	return r
end


# loops over matrix rows
function evalTensor4(mat1::Array{Float64,2},mat2::Array{Float64,2},mat3::Array{Float64,2},mat4::Array{Float64,2},c::Vector)

	n1,m1 = size(mat1)
	n2,m2 = size(mat2)
	n3,m3 = size(mat3)
	n4,m4 = size(mat4)

	row_offset1 = 0
	col_offset1 = 0
	factor1     = 0.0
	row_offset2 = 0
	col_offset2 = 0
	factor2     = 0.0
	row_offset3 = 0
	col_offset3 = 0
	factor3     = 0.0

	r = zeros(n1 * n2 * n3 * n4)

	# loop over rows of 1
	for row_idx1 in 1:n1
		row_offset1 = (row_idx1-1) * n2

		# loop over rows of 2
		for row_idx2 in 1:n2
			row_offset2 = (row_offset1 + (row_idx2-1)) * n3

			# loop rows 3
			for row_idx3 in 1:n3
				row_offset3 = (row_offset2 + (row_idx3-1)) * n4
				# println(row_offset2 + row_idx3)

				for row_idx4 in 1:n4

					# loop over cols of 1
					for col_idx1 in 1:m1
						col_offset1 = (col_idx1-1) * m2
						factor1 = mat1[row_idx1,col_idx1]

						# loop over cols of 2
						for col_idx2 in 1:m2
							col_offset2 = (col_offset1 + (col_idx2-1)) * m3
							factor2 = factor1 * mat2[row_idx2,col_idx2]

							# cols mat3
							for col_idx3 in 1:m3
								col_offset3 = (col_offset2 + (col_idx3-1)) * m4
								factor3 = factor2 * mat3[row_idx3,col_idx3]

								for col_idx4 in 1:m4

									r[row_offset3 + row_idx4] += factor3 * mat4[row_idx4,col_idx4] * c[col_offset3 + col_idx4]
								end
							end
						end
					end
				end
			end
		end
	end
	return r
end

function evalTensor4(mat1::SparseMatrixCSC{T, Int64}, mat2::SparseMatrixCSC{T, Int64}, mat3::SparseMatrixCSC{T, Int64}, mat4::SparseMatrixCSC{T, Int64}, c::Vector{T}) where T

	if size(mat1)[2] > 1
		error("only doing column vector so far, sorry")
	end

	m1 = length(mat1)
	m2 = length(mat2)
	m3 = length(mat3)
	m4 = length(mat4)
	offset1= 0
	factor1= 0.0
	offset2= 0
	factor2= 0.0
	offset3= 0
	factor3= 0.0

	r = 0.0

	# loop over non-zeros in mat1
	for idx1 in mat1.rowval
		offset1 = (idx1-1) * m2
		factor1 = mat1[idx1]

		# loop over non-zeros in mat2
		for idx2 in mat2.rowval
			offset2 = (offset1 + (idx2-1)) * m3
			factor2 = factor1 * mat2[idx2]

			# cols mat3
			for idx3 in mat3.rowval
				offset3 = (offset2 + (idx3-1)) * m4
				factor3 = factor2 * mat3[idx3]

				for idx4 in mat4.rowval
					r += factor3 * mat4[idx4] * c[offset3 + idx4]
				end
			end
		end
	end
	return r
end



function evalTensor(bs::Dict{Int64, SparseMatrixCSC{T, Int64}}, c::Array{T}) where T

	n = length(bs)
	if n == 2
		evalTensor2(bs[2],bs[1],c)
	elseif n==3
		evalTensor3(bs[3],bs[2],bs[1],c)
	elseif n==4
		evalTensor4(bs[4],bs[3],bs[2],bs[1],c)
	else 
		throw(ErrorException("only dims 2,3,4 implemented so far"))
	end
end






