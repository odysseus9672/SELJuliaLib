
"""
	hyperBess_i(d, l, x)

Compute the modified hyperspherical bessel function of the first kind.
Symbolically, this is equal to:
	`besseli(l + d/2 - 1, x) / x^(d/2 - 1)`
with the removable singularity at `x=0` filled in for `l >= 0`.
`hyperBess_i(d, l, k*x)` is an eigenfunction of the `d`-dimensional
Laplacian operator with eigenvalue `-k^2` that corresponds to the
hyperspherical harmonic with eigenvalue `-l*(l+d-2)`.
"""
#Float64 optimized implementation
function hyperBess_i(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_i: only positive values of d are supported")
	end

	if l < 0.0
		error("hyperBess_i: only positive values of l are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > abs(d) * 1e-12
		result = besseli(baseorderplusl, x)
		result *= x^(-baseorder)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = muladd(0.25 * x, x / (baseorderplusl + 1.0), 1.0)
		result *= (0.5*x)^l / ( 2.0^baseorder * gamma(baseorderplusl + 1.0) )
	else
		error("hyperBess_i: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_i

#Generic implementation
function hyperBess_i(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_i: only positive values of d are supported")
	end

	if l < 0
		error("hyperBess_i: only positive values of l are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = l + baseorder

	if abs(x) > abs(d) * 1e-12
		result = besseli(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder)

	else #Do a first order approximation based on the definition
		result = 1 + x * x / (4 * (baseorderplusl + 1))
		result *= (x/2)^l / (2^baseorder * gamma(baseorderplusl + 1))
	end # if abs(x) > abs(d) / 1e12

	return result
end #hyperBess_i


"""
	hyperBess_iprime(d, l, x)

Compute the derivative of the modified hyperspherical Bessel function
of the first kind for
with respect to `x`. Symbolically, this is equal to:
	`hyperBess_i(d, l - 1, x) / 2
	 +hyperBess_i(d, l + 1, x) / 2
	 -(d/2 - 1) * hyperBess_i(d, l, x) / x`
"""
#Float64 optimized implementation
function hyperBess_iprime(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_iprime: only positive values of d are supported")
	end

	if l < 0.0
		error("hyperBess_iprime: only positive values of l are supported")
	end

	if abs(l) < eps(Float64) #send to i0 for processing
		return hyperBess_i0prime( d, x )
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > abs(d) / 1e12
		result = 0.5 * x * (besseli(baseorderplusl - 1.0, x)
			+ besseli(baseorderplusl + 1.0, x))
		result = muladd(-baseorder, besseli(baseorderplusl, x), result)
		result *= x^(-baseorder)

	elseif 0.0 <= x #Do a three term approximation based on the definition
		xsqr::Float64 = x*x
		result = muladd(0.25 * (l+4.0)/(baseorderplusl + 2.0), xsqr, (l+2) )
		result /= 4.0 * (baseorderplusl + 1.0)
		result = muladd(result, xsqr, l)
		result *= z^(l-1) / (2.0^baseorderplusl * gamma(baseorderplusl + 1.0))

	else
		error("hyperBess_i0prime: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_i0prime


#Generic implementation
function hyperBess_iprime(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_i0prime: only positive values of d are supported")
	end

	if l < 0
		error("hyperBess_iprime: only positive values of l are supported")
	end

	if abs(l) < eps(Float64) #send to i0 for processing
		return hyperBess_i0prime( d, x )
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = convert(Tp, l) + baseorder

	if abs(x) > abs(d) / 1e12
		result = x * (besseli(real(baseorderplusl) - 1, x)
			+ besseli(real(baseorderplusl) + 1, x)) / 2
		result -= baseorder * besseli(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder) * gamma(baseorder + 1)

	else
		xsqr = x * x
		result = xsqr * (l + 4)/(baseorderplusl + 2) + (l+2)
		result /= 4 * (baseorderplusl + 1)
		result = xsqr * result + l
		result *= z^(l-1) / (2^baseorderplusl * gamma(baseorderplusl + 1))
	end # if x > d / 1e12

	return result
end #hyperBess_i0prime

