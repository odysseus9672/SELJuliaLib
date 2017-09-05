
"""
    hyperBess_y0(d, x)

Compute the hyperspherical bessel function of the second kind for symmetrical states.
Symbolically, this is equal to:
	`bessely(d/2 - 1, x) / x^(d/2-1)`
with the removable singularity at `x=0` filled in. `hyperBess_y0(d, k*x)` is
a hyperspherically symmetric eigenfunction of the `d`-dimensional Laplacian operator
with eigenvalue `-k^2`. For `d = 1`, this function is identical to `sin(x)`
multiplied by `sqrt(2/pi)`, for `d=2` it corresponds to `bessely(0,x)`,
and for `d=3` it gives the spherical Bessel functions of the second kind
multiplied by `sqrt(2/pi)`.
"""
#Float64 optimized implementation
function hyperBess_y0(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_y0: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > 0.0
		result = bessely(baseorder, x) * x^(-baseorder)

	else
		error("hyperBess_y0: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_y0


#Generic implementation
function hyperBess_y0(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_y0: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > 0
		result = bessely(real(baseorder), x) * x^(-baseorder)

	else
		error("hyperBess_y0: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_y0


"""
	hyperBess_y0prime(d, x)

Compute the derivative of the hyperspherical Bessel function of the second kind for
symmetrical states with respect to `x`. Symbolically, this is equal to:
	`-hyperBess_y(d, 1, x)`
according to DLMF equation 10.6.6 http://dlmf.nist.gov/10.6#E6
"""
#Float64 optimized implementation
function hyperBess_y0prime(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_y0prime: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > 0.0
		result = - bessely(baseorder + 1.0, x) * x^(-baseorder)

	else
		error("hyperBess_y0prime: non-positive x outside of domain.")
	end # if x > 0.0

	return result
end#::Float64 #hyperBess_y0prime


#Generic implementation
function hyperBess_y0prime(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_y0prime: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > 0
		result = - bessely(real(baseorder) + 1, x) * x^(-baseorder)

	else
		error("hyperBess_y0prime: x = 0 outside of domain.")
	end # if abs(x) > 0

	return result
end #hyperBess_y0prime


"""
	hyperBess_y0pprime(d, x)

Compute the 2nd derivative of the hyperspherical Bessel function of the second kind for
symmetrical states with respect to `x`. Symbolically, this is equal to:
	`- hyperBess_y(d, 0, x) / 2
		+ hyperBess_y(d, 2, x) / 2
		+ (d/2 - 1) * hyperBess_y(d, 1, x) / x`
This function is reserved for internal use for finding zeros of `hyperBess_y0prime`.
"""
#Float64 optimized implementation
function hyperBess_y0pprime(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_y0pprime: only positive values of d are supported")
	end

	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x >= 0.0
		result = 0.5*x*(bessely(baseorder+2.0, x) - bessely(baseorder, x))
		result = muladd(baseorder, bessely(baseorder + 1.0, x), result)
		result *= x^(baseorder + 1.0)

	else
		error("hyperBess_y0pprime: non-positive x outside of domain.")
	end

	return result
end#::Float64 #hyperBess_y0pprime


#Generic implementation
function hyperBess_y0pprime(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_y0pprime: only positive values of d are supported")
	end
	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > 0
		result = x*(bessely(real(baseorder)+2, x) - bessely(real(baseorder), x)) / 2
		result += baseorder * bessely(real(baseorder) + 1, x)
		result *= 2^baseorder * gamma(baseorder + 1) * x^(baseorder + 1)

	else
		error("hyperBess_y0pprime: x = 0 outside of domain.")
	end #if abs(x) > 0

	return result
end #hyperBess_y0pprime


"""
	hyperBess_y0zeros(d, mmax)

Compute the first mmax zeros of the hyperspherical Bessel function of the second kind for
symmetrical states.
"""
#Float64 optimized version
function hyperBess_y0zeros( d::Float64, mmax::Int )
	if mmax <= 0
		error("hyperBess_y0zeros: mmax must be positive")
	end

	if d <= 0.0
		error("hyperBess_y0zeros: d must be positive")
	end

	baseorder::Float64 = muladd(0.5, d, -1.0)
	#empirical function that gets near the first zero
	guess::Float64 = sqrt((baseorder + 2.0)^2 - 1.0)

	result::Array{Float64,1} = Array{Float64}(mmax)
	converged::Bool = true

	function f(x)
		return hyperBess_y0(d, x)
	end

	function fp(x)
		return hyperBess_y0prime(d, x)
	end

	for m in 1:mmax
		result[m], converged = findrootNewtonF64toF64(f, fp, guess)
		if !converged
			error("hyperBess_y0zeros: failed to converge on zero $m with d=$d")
		end

		guess = result[m] + pi
	end

	return result
end#::Array{Float64,1} #hyperBess_y0zeros


#Generic implementation
function hyperBess_y0zeros( d::Real, mmax::Int, exampleval::Tp ) where Tp <: Number
	if mmax <= 0
		error("hyperBess_y0zeros: mmax must be positive")
	end

	if d <= 0
		error("hyperBess_y0zeros: d must be positive")
	end

	baseorder = convert(Tp, d)/2 - 1
	#empirical function that gets near the first zero
	guess = sqrt((baseorder + 2)^2 - 1)

	result = Array{Tp}(mmax)
	converged::Bool = true

	function f(x::Tp)
		hyperBess_y0(d, x)
	end

	function fp(x::Tp)
		hyperBess_y0prime(d, x)
	end

	for m in 1:mmax
		result[m], converged = findrootNewton(f, fp, guess)
		if !converged
			error("hyperBess_y0zeros: failed to converge on zero $m with d=$d")
		end

		guess += pi
	end

	return result
end #hyperBess_y0zeros