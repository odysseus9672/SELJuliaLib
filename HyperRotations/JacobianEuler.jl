"""
	jacobianSONEuler(N, thetas)

	The space of rotation matrices represents SO(N), and the geometry of
	SO(N) is isomorphic to the product of the surfaces of n-dimensional spheres,
	with n going from 2 (a circle) to N. Integrating the Jacobian over the
	canonical values that Euler angles can take yields the product of the areas
	of those spheres. In other words, it relates the coordinate valome in angles
	to the true volume of the space. Another name for the true volume, as its used
	here, is the Haar measure. It serves as the density function for rotation
	matrices randomly selected with uniform probability.
"""

function jacobianSONEuler(N::Int, thetas::Vector{Tp}) where Tp<:Real
	lenAngs = length(thetas)

	if N * (N - 1) != 2 * lenAngs
		error("The number of dimensions (N) must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	if lenAngs < 1
		return NaN
	end

	Jacobian::Tp = 1
	lenUvecAngs::Int = 2
	ith::Int = 2 #begin with the second vector
	while ith <= lenAngs
		vecEnd = ith + lenUvecAngs - 1
		ith += 1 #The first angle in the unit vector contributes a factor of 1

		curPow::Int = 1
		while ith <= vecEnd
			Jacobian *= abs(sin(thetas[ith]))^curPow

			curPow += 1
			ith += 1
		end

		lenUvecAngs += 1
	end

	return abs(Jacobian)
end
