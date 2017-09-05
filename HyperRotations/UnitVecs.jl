"""
	unitVecFromPolarAngles( thetas )

	Returns the unit vector defined by the array of angles, `thetas`, that
	has length one greater than `thetas`. N-dimensional polar
	unit vectors are defined recursively by the relation:

	rhat_N = sin(theta_{N-1}) * rhat_{N-1} + cos(theta_{N-1}) * ehat_{N}

	where rhat_N is the general N-dimensional unit vector, theta_{N-1} is
	the last polar angle in a list of length N-1 that defines rhat_N,
	and ehat_{N} is the Nth orthonormal basis vector. In 2-d this gives:

	rhat_2 = [ sin(theta), cos(theta) ].

	In 3-d it's:

	rhat_3 = [ sin(phi) * sin(theta), cos(phi) * sin(theta), cos(theta) ].

	Referring to an angle as polar, then, means that it resembles a rotation
	from the higher index basis vector into the lower one, just like the polar
	angle in the classical parametrization of 3-d vectors."""

function unitVecFromPolarAngles( thetas::Vector{Tp} ) where Tp<:Real
	numthetas::Int = length(thetas)
	uvec = ones(Tp, numthetas + 1)

	i = 1
	while i <= numthetas
		s, c = sincos(thetas[i])
		uvec[1:i] *= s
		i += 1
		uvec[i] *= c
	end

	return uvec
end


"""normalizeVec!( v )
	Divides the argument, in place, by its Euclidean norm
	to give it length 1."""

function normalizeVec!( v::Vector{Tp} ) where Tp<:Real
	v[1:end] = v ./ HypotVec(v)
	return
end


"""
	normalizeVec( v )

	normalizeVec returns a copy of the argument divided by its Euclidean
	norm, giving a unit vector."""

function normalizeVec( v::Vector{Tp} ) where Tp<:Real
	return v ./ hypotVec(v)
end


"""
	polarAnglesFromUnitVec( uvec, zeroTol = 8 * eps() )

	Inverts unitVecFromPolarAngles, returning the angles in the pure polar
	parametrization of its argument. The argument `zeroTol` specifies magnitude
	below which values will be treated as exactly zero. This is done to avoid
	divide by zero errors and reduce numerical instability.
"""

function polarAnglesFromUnitVec( uvec::Vector{Tp}, zeroTol::Tp = 8 * eps(Tp) ) where Tp<:Real
	veclen = length(uvec)

	if veclen >= 2
		numthetas = veclen - 1
		thetas = Array{Tp}(numthetas)

		#Handle the purely polar angles
		ith = numthetas
		denom::Tp = 1.0
		iv = veclen
		while iv > 1 && denom > zeroTol
			th = acos(uvec[iv])
			thetas[ith] = th
			denom *= sin(th)

			iv -= 1
			ith -= 1
		end

		#Set any remaining purely polar angles to zero
		while ith > 1
			thetas[ith] = 0.0
			ith -= 1
		end

		#Set the azimuthal angle - Note that this is intentionally cos = y, sin = x
		thetas[1] = atan2(uvec[1], uvec[2])

	else
		thetas = Array{Tp}(0)
	end

	return thetas
end
