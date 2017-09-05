
"""
	hyperBess_j(d, l, x)

Compute the hyperspherical bessel function of the first kind.
Symbolically, this is equal to:
	`besselj(l + d/2 - 1, x) / x^(d/2 - 1)`
with the removable singularity at `x=0` filled in for `l >= 0`.
`hyperBess_j(d, l, k*x)` is an eigenfunction of the `d`-dimensional
Laplacian operator with eigenvalue `-k^2` that corresponds to the
hyperspherical harmonic with eigenvalue `-l*(l+d-2)`. For `d=2` this
function reproduces `besselj(l,x)`, and for `d=3` it gives the
spherical Bessel functions of the first kind multiplied by `sqrt(2/pi)`.
"""
#Float64 optimized implementation
function hyperBess_j(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_j: only positive values of d are supported")
	end

	if l < 0.0
		error("hyperBess_j: only positive values of l are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > abs(d) * 1e-12
		result = besselj(baseorderplusl, x)
		result *= x^(-baseorder)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = muladd(-0.25 * x, x / (baseorderplusl + 1.0), 1.0)
		result *= (0.5*x)^l / ( 2.0^baseorder * gamma(baseorderplusl + 1.0) )
	else
		error("hyperBess_j: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_j

#Generic implementation
function hyperBess_j(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_j: only positive values of d are supported")
	end

	if l < 0
		error("hyperBess_j: only positive values of l are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = l + baseorder

	if abs(x) > abs(d) * 1e-12
		result = besselj(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder)

	else #Do a first order approximation based on the definition
		result = 1 - x * x / (4 * (baseorderplusl + 1))
		result *= (x/2)^l / (2^baseorder * gamma(baseorderplusl + 1))
	end # if abs(x) > abs(d) / 1e12

	return result
end #hyperBess_j


"""
	hyperBess_jprime(d, l, x)

Compute the derivative of the hyperspherical Bessel function of the first kind for
with respect to `x`. Symbolically, this is equal to:
	`hyperBess_j(d, l - 1, x) / 2
	 -hyperBess_j(d, l + 1, x) / 2
	 -(d/2 - 1) * hyperBess_j(d, l, x) / x`
"""
#Float64 optimized implementation
function hyperBess_jprime(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_jprime: only positive values of d are supported")
	end

	if l < 0.0
		error("hyperBess_jprime: only positive values of l are supported")
	end

	if abs(l) < eps(Float64) #send to j0 for processing
		return hyperBess_j0prime( d, x )
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > abs(d) / 1e12
		result = 0.5 * x * (besselj(baseorderplusl - 1.0, x)
			- besselj(baseorderplusl + 1.0, x))
		result = muladd(-baseorder, besselj(baseorderplusl, x), result)
		result *= x^(-baseorder)

	elseif 0.0 <= x #Do a three term approximation based on the definition
		xsqr::Float64 = x*x
		result = muladd(0.25 * (l+4.0)/(baseorderplusl + 2.0), xsqr, -(l+2) )
		result /= 4.0 * (baseorderplusl + 1.0)
		result = muladd(result, xsqr, l)
		result *= z^(l-1) / (2.0^baseorderplusl * gamma(baseorderplusl + 1.0))

	else
		error("hyperBess_j0prime: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_j0prime


#Generic implementation
function hyperBess_jprime(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_j0prime: only positive values of d are supported")
	end

	if l < 0
		error("hyperBess_jprime: only positive values of l are supported")
	end

	if abs(l) < eps(Float64) #send to j0 for processing
		return hyperBess_j0prime( d, x )
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = convert(Tp, l) + baseorder

	if abs(x) > abs(d) / 1e12
		result = x * (besselj(real(baseorderplusl) - 1, x)
			- besselj(real(baseorderplusl) + 1, x)) / 2
		result -= baseorder * besselj(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder) * gamma(baseorder + 1)

	else
		xsqr = x * x
		result = xsqr * (l + 4)/(baseorderplusl + 2) - (l+2)
		result /= 4 * (baseorderplusl + 1)
		result = xsqr * result + l
		result *= z^(l-1) / (2^baseorderplusl * gamma(baseorderplusl + 1))
	end # if x > d / 1e12

	return result
end #hyperBess_j0prime

