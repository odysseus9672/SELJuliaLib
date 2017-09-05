
"""
    hyperBess_i0(d, x)

Compute the modified hyperspherical bessel function of the first kind for
symmetrical states.
Symbolically, this is equal to:
	`besseli(d/2 - 1, x) / x^(d/2-1)`
with the removable singularity at `x=0` filled in. `hyperBess_i0(d, k*x)` is
a hyperspherically symmetric eigenfunction of the `d`-dimensional Laplacian operator
with eigenvalue `k^2`.
"""
#Float64 optimized implementation
function hyperBess_i0(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_i0: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > abs(d) * 1e-12
		result = besseli(baseorder, x) * x^(-baseorder)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = muladd(0.25 * x / (baseorder + 1.0), x, 1.0)
		result /= 2.0^baseorder * gamma(0.5*d)
	else
		error("hyperBess_i0: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_i0


#Generic implementation
function hyperBess_i0(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_i0: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > abs(d) * 1e-12
		result = besseli(real(baseorder), x) * x^(-baseorder)

	else #Do a first order approximation based on the definition
		result = 1 + x * x / (4 * (baseorder + 1))
		result /= 2^baseorder * gamma(baseorder + 1)
	end # if abs(x) > abs(d) / 1e12

	return result
end #hyperBess_i0


"""
	hyperBess_i0prime(d, x)

Compute the derivative of the modified hyperspherical Bessel function
of the first kind for
symmetrical states with respect to `x`. Symbolically, this is equal to:
	`hyperBess_i(d, 1, x)`
according to DLMF equation 10.29.4 http://dlmf.nist.gov/10.29#E4
"""
#Float64 optimized implementation
function hyperBess_i0prime(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_i0prime: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > abs(d) / 1e12
		result = besseli(baseorder + 1.0, x) * x^(-baseorder)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = x * muladd(0.5 * x, x/(baseorder + 2.0), 1.0) / d
		result /= 2.0^baseorder * gamma(baseorder + 1.0)

	else
		error("hyperBess_i0prime: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_i0prime


#Generic implementation
function hyperBess_i0prime(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_i0prime: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if x > abs(d) / 1e12
		result = besseli(real(baseorder) + 1, x) * x^(-baseorder)

	else
		result = x * ( x*x/(2*(baseorder+2)) + 1 ) / d
		result /= 2^baseorder * gamma(baseorder + 1)
	end # if x > d / 1e12

	return result
end #hyperBess_i0prime