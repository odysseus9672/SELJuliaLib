"""
	matrixSONLie( N, thetas )

	matrixSONLie returns an SO(N) rotation matrix that comes from constructing
	an antisymmetric matrix, called a generator of a rotation in terms of Lie algebra
	(group theory) , from the angles in thetas. The rotation matrix is the matrix
	exponential of the generator, R = expm(G). The convention for forming the generator
	is chosen so that if all thetas are positive the diagonal above the main is negative,
	and the diagonals alternate signs going to the upper right corner. For example,
	the generator for SO(3) is formed as:
	G = Tp[  0         -thetas[1]  thetas[2];
		     thetas[1]  0         -thetas[3];
		    -thetas[2]  thetas[3]  0         ].
	This function is, roughly, the generalization of the axis and angle formalism in
	SO(3). The translation is:
		theta = hypot(thetas)
		axis[1] = thetas[3] / theta
		axis[2] = thetas[2] / theta
		axis[3] = thetas[1] / theta
"""

function matrixSONLie( N::Integer, thetas::Vector{Tp} ) where Tp<:Real
	lenAngs = length(thetas)

	if N * (N - 1) != 2 * lenAngs
		error("The number of dimensions (N) must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	if lenAngs < 1
		return Array{Tp}((0,0))
	end

	generator = Array{Tp}(N, N)
	k = 1
	for fromidx in 1:N
		generator[fromidx,fromidx] = zero(Tp)
		sgn = -one(Tp)
		for toidx in (fromidx+1):N
			val = sgn * thetas[k]
			generator[i, j] = val
			generator[j, i] = -val
			k += 1
			sgn *= -one(Tp)
		end
	end

	return expm(generator)
end
