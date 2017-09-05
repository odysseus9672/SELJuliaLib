
"""
    hyperBess_k0(d, x)

Compute the modified hyperspherical bessel function of the second kind for
symmetrical states.
Symbolically, this is equal to:
	`besselk(d/2 - 1, x) / x^(d/2-1)`
with the removable singularity at `x=0` filled in. `hyperBess_k0(d, k*x)` is
a hyperspherically symmetric eigenfunction of the `d`-dimensional Laplacian operator
with eigenvalue `k^2`.
"""
#Float64 optimized implementation
function hyperBess_k0(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_k0: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > 0.0
		result = besselk(baseorder, x) * x^(-baseorder)

	else
		error("hyperBess_k0: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_k0


#Generic implementation
function hyperBess_k0(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_k0: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > 0
		result = besselk(real(baseorder), x) * x^(-baseorder)

	else
		error("hyperBess_k0: x = 0 outside of domain."))
	end # if abs(x) > 0

	return result
end #hyperBess_k0


"""
	hyperBess_k0prime(d, x)

Compute the derivative of the modified hyperspherical Bessel function
of the first kind for
symmetrical states with respect to `x`. Symbolically, this is equal to:
	`hyperBess_k(d, 1, x)`
according to DLMF equation 10.29.4 http://dlmf.nist.gov/10.29#E4
"""
#Float64 optimized implementation
function hyperBess_k0prime(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_k0prime: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > 0.0
		result = besselk(baseorder + 1.0, x) * x^(-baseorder)

	else
		error("hyperBess_k0prime: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_k0prime


#Generic implementation
function hyperBess_k0prime(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_k0prime: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > 0
		result = besselk(real(baseorder) + 1, x) * x^(-baseorder)

	else
		error("hyperBess_k0prime: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_k0prime