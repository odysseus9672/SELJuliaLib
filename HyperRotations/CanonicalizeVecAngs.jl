
function canonicalizeUnitVecAngles!(thetas::Vector{Tp}) where Tp<:Real
	#Bring all angles into (-pi, pi]
	thetas[1:end] = map(mod2pi, thetas)
	thetas[1:end] -= 2*pi .* ceil(thetas ./ (2*pi) .- 0.5)

	#Bring the polar angles into their canonical range: [0, pi]
	ith = length(thetas)
	while ith > 1
		if thetas[ith] >= 0.0 #already canonical
			ith -= 1
			continue
		end

		#Need to adjust theta[ith], and all lower angles
		thetas[ith] = -thetas[ith]
		jth = ith - 1
		while jth > 1
			#Shift lower angles without leaving (-pi, pi]
			if thetas[jth] > 0.0
				thetas[jth] = pi - thetas[jth]
			else
				thetas[jth] = -pi - thetas[jth]
			end

			jth -= 1
		end

		#Shift the azumuthal angle, again without leaving (-pi, pi]
		if thetas[1] > 0.0
			thetas[1] = thetas[1] - pi
		else
			thetas[1] = thetas[1] + pi
		end

		ith -= 1
	end

	#Shift the azumthal ingle into [0.0, 2*pi)
	thetas[1] = mod2pi(thetas[1])
	return
end



"""
	canonicalizeUnitVecAngles(thetas)

	canonicalizeUnitVecAngles brings the angles that specify a unit vector into
	their canonical ranges without altering the components of the unit vector the
	angles represent. The first angle is an 'azimuthal' polar angle and has
	the range [0, 2*pi). All subsequent angles are purely polar, and have the
	range [0, pi].
"""

function canonicalizeUnitVecAngles(thetas::Vector{Tp}) where Tp<:Real
	Cthetas = copy(thetas)
	canonicalizeUnitVecAngles!(Cthetas)

	return cThetas
end