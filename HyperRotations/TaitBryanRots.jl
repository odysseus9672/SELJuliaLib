
function rotateSONTaitBryan!(v::Vector{Tp}, thetas::Vector{Tp}) where Tp<:Real
	lenAngs = length(thetas)

	if N * (N - 1) != 2 * lenAngs
		error("The number of dimensions (N) must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	if lenAngs < 1
		return Array{Tp}((0,0))
	end

	coordPairs = [ (1, j) for j in 2:N ]
	for i in 2:(N-1)
		coordPairs = vcat(coordPairs, [ (i, j) for j in (i+1):N ])
	end

	rotateSeq!(v, thetas, coordPairs)

	return
end


"""
	rotateSONTaitBryan(v, thetas)

	Rotates the vector it is given in the same way as the
	matrix that corresponds to the angles in thetas. The pure polar representation
	of rotation matrices is defined in the documentation for
	matrixSONEuler.
"""

function rotateSONTaitBryan(v::Vector{Tp}, thetas::Vector{Tp}) where Tp<:Real
	result = copy(v)
	N = length(v)
	if N * (N - 1) != 2 * length(thetas)
		error("The number of dimensions (N) of v must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	coordPairs = [ (1, j) for j in 2:N ]
	for i in 2:(N-1)
		coordPairs = vcat(coordPairs, [ (i, j) for j in (i+1):N ])
	end

	rotateSeq!(result, thetas, coordPairs)
	return result
end
