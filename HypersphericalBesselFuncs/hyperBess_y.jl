
"""
	hyperBess_y(d, l, x)

Compute the hyperspherical bessel function of the second kind.
Symbolically, this is equal to:
	`bessely(l + d/2 - 1, x) * 2^(d/2 - 1)`
with the removable singularity at `x=0` filled in for `l >= 0`.
`hyperBess_y(d, l, k*x)` is an eigenfunction of the `d`-dimensional
Laplacian operator with eigenvalue `-k^2` that corresponds to the
hyperspherical harmonic with eigenvalue `-l*(l+d-2)`. For `d=2` this
function reproduces `bessely(l,x)`, and for `d=3` it gives the
spherical Bessel functions of the second kind multiplied by `sqrt(2/pi)`.
"""
#Float64 optimized implementation
function hyperBess_y(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_y: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > 0.0
		result = bessely(baseorderplusl, x)
		result *= x^(-baseorder)

	else
		error("hyperBess_y: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_y

#Generic implementation
function hyperBess_y(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_y: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = l + baseorder

	if abs(x) > 0
		result = bessely(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder)

	else
		error("hyperBess_y: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_y


"""
	hyperBess_yprime(d, l, x)

Compute the derivative of the hyperspherical Bessel function of the first kind for
with respect to `x`. Symbolically, this is equal to:
	`hyperBess_y(d, l - 1, x) / 2
	 -hyperBess_y(d, l + 1, x) / 2
	 -(d/2 - 1) * hyperBess_y(d, l, x) / x`
"""
#Float64 optimized implementation
function hyperBess_yprime(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_yprime: only positive values of d are supported")
	end

	if abs(l) < eps(Float64) #send to j0 for processing
		return hyperBess_y0prime( d, x )
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > 0.0
		result = 0.5 * x * (bessely(baseorderplusl - 1.0, x)
			- bessely(baseorderplusl + 1.0, x))
		result = muladd(-baseorder, bessely(baseorderplusl, x), result)
		result *= x^(-baseorder)

	else
		error("hyperBess_y0prime: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_y0prime


#Generic implementation
function hyperBess_yprime(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_y0prime: only positive values of d are supported")
	end

	if abs(l) < eps(Float64) #send to j0 for processing
		return hyperBess_y0prime( d, x )
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = convert(Tp, l) + baseorder

	if abs(x) > 0
		result = x * (bessely(real(baseorderplusl) - 1, x)
			- bessely(real(baseorderplusl) + 1, x)) / 2
		result -= baseorder * bessely(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder) * gamma(baseorder + 1)

	else
		error("hyperBess_yprime: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_y0prime

