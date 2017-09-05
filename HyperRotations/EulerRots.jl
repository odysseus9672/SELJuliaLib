
function rotateSONEuler!(v::Vector{Tp}, thetas::Vector{Tp}) where Tp<:Real
	N = length(v)
	if N * (N - 1) != 2 * length(thetas)
		error("The number of dimensions (N) of v must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	coordPairs = Array{Tp}((length(thetas), 2) )
	i = 1
	for maxcoord in 1:N
		for fromcoord in 1:(maxcoord - 1)
			coordPairs[i,1:2] = [fromcoord, fromcoord + 1]
			i += 1
		end
	end

	rotateSeq!(v, thetas, coordPairs)
end


"""
	rotateSONEuler(v, thetas)

	rotateSONEuler rotates the vector it is given in the same way as the
	matrix that corresponds to the angles in thetas. The pure polar representation
	of rotation matrices is defined in the documentation for
	matrixSONEuler.
"""

function rotateSONEuler(v::Vector{Tp}, thetas::Vector{Tp}) where Tp<:Real
	result = copy(v)
	N = length(v)
	if N * (N - 1) != 2 * length(thetas)
		error("The number of dimensions (N) of v must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	coordPairs = Array{Tp}((length(thetas), 2) )
	i = 1
	for maxcoord in 1:N
		for fromcoord in 1:(maxcoord - 1)
			coordPairs[i,1:2] = [fromcoord, fromcoord + 1]
			i += 1
		end
	end

	rotateSeq!(result, thetas, coordPairs)
	return result
end
