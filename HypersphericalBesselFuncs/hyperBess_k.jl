
"""
	hyperBess_k(d, l, x)

Compute the modified hyperspherical bessel function of the second kind.
Symbolically, this is equal to:
	`besselk(l + d/2 - 1, x) / x^(d/2 - 1)`
with the removable singularity at `x=0` filled in for `l >= 0`.
`hyperBess_k(d, l, k*x)` is an eigenfunction of the `d`-dimensional
Laplacian operator with eigenvalue `-k^2` that corresponds to the
hyperspherical harmonic with eigenvalue `-l*(l+d-2)`.
"""
#Float64 optimized implementation
function hyperBess_k(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_k: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > 0.0
		result = besselk(baseorderplusl, x)
		result *= x^(-baseorder)

	else
		error("hyperBess_k: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_k

#Generic implementation
function hyperBess_k(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_k: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = l + baseorder

	if abs(x) > 0
		result = besselk(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder)

	else
		error("hyperBess_k: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_k


"""
	hyperBess_kprime(d, l, x)

Compute the derivative of the modified hyperspherical Bessel function
of the second kind for
with respect to `x`. Symbolically, this is equal to:
	`hyperBess_k(d, l - 1, x) / 2
	 +hyperBess_k(d, l + 1, x) / 2
	 -(d/2 - 1) * hyperBess_k(d, l, x) / x`
"""
#Float64 optimized implementation
function hyperBess_kprime(d::Float64, l::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_kprime: only positive values of d are supported")
	end

	if abs(l) < eps(Float64) #send to k0 for processing
		return hyperBess_k0prime( d, x )
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplusl::Float64 = l + baseorder

	if x > 0.0
		result = 0.5 * x * (besselk(baseorderplusl - 1.0, x)
			+ besselk(baseorderplusl + 1.0, x))
		result = muladd(-baseorder, besselk(baseorderplusl, x), result)
		result *= x^(-baseorder)

	else
		error("hyperBess_k0prime: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_k0prime


#Generic implementation
function hyperBess_kprime(d::Real, l::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_k0prime: only positive values of d are supported")
	end

	if l < 0
		error("hyperBess_kprime: only positive values of l are supported")
	end

	if abs(l) < eps(Float64) #send to k0 for processing
		return hyperBess_k0prime( d, x )
	end

	baseorder = convert(Tp, d) / 2 - 1
	baseorderplusl = convert(Tp, l) + baseorder

	if abs(x) > 0
		result = x * (besselk(real(baseorderplusl) - 1, x)
			+ besselk(real(baseorderplusl) + 1, x)) / 2
		result -= baseorder * besselk(real(baseorderplusl), x)
		result *= (x/2)^(-baseorder) * gamma(baseorder + 1)

	else
		error("hyperBess_kprime: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_k0prime

