"""
	versin(t)

	Calculates the versed sine of `t`.
	The definition of versin is:
	``\\mathrm{versin}(x) = 1 - \\cos(x).``
	By trig identies, it can be shown that:
	``\\mathrm{versin}(x) = 2  sin^2\\left(\\frac{x }{ 2}\\right)^2.``
	This function uses the second version when `x = mod2pi(t) <= 1`,
	the first otherwise.
"""

function versin(t::Tp) where Tp<:Real
	x = mod2pi(t)
	return x <= one(Tp) ? 2 * sin(x/2)^2 : one(Tp) - cos(x)
end



"""
	coversin(t)
	Calculates the coversed sine (cvs) of its argument.
	The definition of cvs is:
	``\\mathrm{coversin}(x) = 1 - \\sin(x).``
	By trig identies, it can be shown that:
	``\\mathrm{coversin}(x) = 2 \\mathrm{sin}^2\\left(\\frac{2x - \\pi}{4}\\right)^2.``
	This function uses the second version when `x = mod2pi(t) > 1`,
	the first otherwise.
"""

function coversin(t::Tp) where Tp<:Real
	x = mod2pi(t)
	return x <= one(Tp) ? one(Tp) - sin(x) : 2 * sin(x/2 - pi/4)^2
end
