
function canonicalizeSONEuler!(N::Integer, thetas::Vector{Tp}) where Tp<:Real
	lenAngs = length(thetas)

	if N * (N - 1) != 2 * lenAngs
		error("The number of dimensions (N) must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	if lenAngs < 1
		return
	end

	#Bring all angles into (-pi, pi]
	thetas[1:end] = map(mod2pi, thetas)
	thetas[1:end] -= 2*pi .* ceil(thetas ./ (2*pi) .- one(Tp) / 2)

	#initialize the entries of the diagonal sign matrix
	signMatDiag = ones(Tp, N)

	#Begin the main loop
	vecEnd::Int = lenAngs
	vecLen::Int = N - 1
	vecStart = vecEnd - vecLen + 1

	while vecLen > 0
		#=Move the sign matrix from the right of the rotation matrix defined by the
		current unit vector to its left
		=#

		isgn::Int = 1
		jsgn::Int = 2
		iang::Int = vecStart
		while iang <= vecEnd
			if signMatDiag[isgn] == -signMatDiag[jsgn]
				thetas[iang] *= -1
			end

			isgn = jsgn
			jsgn += 1
			iang += 1
		end

		#Cancel the lower right minus sign in signMat, if needed
		iang -= 1
		if signMatDiag[isgn] < 0
			signMatDiag[isgn] = 1
			signMatDiag[1] *= -1

			while iang > 1
				thetas[iang] = pi - thetas[iang]
				iang -= 1
			end

			if thetas[vecStart] > 0
				thetas[vecStart] -= pi
			else
				thetas[vecStart] += pi
			end
		end

		#Canonicalize the polar angles in the current vector
		#iang and isgn are already near the right place due to last loop
		isgn -= 1
		iang = vecEnd
		while iang > vecStart
			if thetas[iang] < 0
				thetas[iang] *= -1
				signMatDiag[isgn] *= -1
				signMatDiag[1] *= -1

				jang::Int = iang - 1
				while jang > vecStart
					thetas[jang] = pi - thetas[jang]

					jang -= 1
				end

				if thetas[vecStart] > 0
					thetas[vecStart] -= pi
				else
					thetas[vecStart] += pi
				end
			end

			iang -= 1
			isgn -= 1
		end

		#Shift the azimuthal angle into [0, 2*pi)
		thetas[vecStart] = mod2pi(thetas[vecStart])

		vecLen -= 1
		vecEnd = vecStart - 1
		vecStart -= vecLen
	end

	return
end


"""
	canonicalizeSONEuler(N, thetas)

	canonicalizeSONEuler brings the angles stored in thetas to all be in
	their canonical ranges. The matrix produced from the angles will be
	identical, up to rounding errors. The angles in thetas are intepreted
	as a list of unit vector angles with dimensions from 2 to N. The
	canonical ranges are defined as the canonical ranges for the individual
	unit vectors. If an N-dimensional unit vector has N-1 angles, theta_i,
	the canonical ranges are:
		theta_1 in [0, 2*pi), and all other theta_i in [0, pi].
	Concretely, for N = 4 the pattern of upper limits is:
	[ 2 * pi,
	  2 * pi, pi,
	  2 * pi, pi, pi ].
"""
function canonicalizeSONEuler(N::Integer, thetas::Vector{Tp}) where Tp<:Real
	result = copy(thetas)

	canonicalizeSONEuler!(N, result)
	return result
end
