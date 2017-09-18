
"""
    hyperBess_j0(d, x)

Compute the hyperspherical bessel function of the first kind for symmetrical states.
Symbolically, this is equal to:
	`besselj(d/2 - 1, x) / x^(d/2-1)`
with the removable singularity at `x=0` filled in. `hyperBess_j0(d, k*x)` is
a hyperspherically symmetric eigenfunction of the `d`-dimensional Laplacian operator
with eigenvalue `-k^2`. For `d = 1`, this function is identical to `cos(x)`
multiplied by `sqrt(2/pi)`, for `d=2` it corresponds to `besselj(0,x)`,
and for `d=3` it gives the spherical Bessel functions of the first kind
multiplied by `sqrt(2/pi)`.
"""
#Float64 optimized implementation
function hyperBess_j0(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_j0: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > abs(d) * 1e-12
		result = besselj(baseorder, x) * x^(-baseorder)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = muladd(-0.25 * x / (baseorder + 1.0), x, 1.0)
		result /= 2.0^baseorder * gamma(0.5*d)
	else
		error("hyperBess_j0: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_j0


#Generic implementation
function hyperBess_j0(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_j0: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if abs(x) > abs(d) * 1e-12
		result = besselj(real(baseorder), x) * x^(-baseorder)

	else #Do a first order approximation based on the definition
		result = 1 - x * x / (4 * (baseorder + 1))
		result /= 2^baseorder * gamma(baseorder + 1)
	end # if abs(x) > abs(d) / 1e12

	return result
end #hyperBess_j0


"""
	hyperBess_j0prime(d, x)

Compute the derivative of the hyperspherical Bessel function of the first kind for
symmetrical states with respect to `x`. Symbolically, this is equal to:
	`-hyperBess_j(d, 1, x)`
according to DLMF equation 10.6.6 http://dlmf.nist.gov/10.6#E6
"""
#Float64 optimized implementation
function hyperBess_j0prime(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_j0prime: only positive values of d are supported")
	end

	result::Float64 = NaN
	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > abs(d) / 1e12
		result = - besselj(baseorder + 1.0, x) * x^(-baseorder)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = x * muladd(0.5 * x, x/(baseorder + 2.0), -1.0) / d
		result /= 2.0^baseorder * gamma(baseorder + 1.0)

	else
		error("hyperBess_j0prime: negative x outside of domain.")
	end # if x > d / 1e12

	return result
end#::Float64 #hyperBess_j0prime


#Generic implementation
function hyperBess_j0prime(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_j0prime: only positive values of d are supported")
	end

	baseorder = convert(Tp, d) / 2 - 1

	if x > abs(d) / 1e12
		result = - besselj(real(baseorder) + 1, x) * x^(-baseorder)

	else
		result = x * ( x*x/(2*(baseorder+2)) - 1 ) / d
		result /= 2^baseorder * gamma(baseorder + 1)
	end # if x > d / 1e12

	return result
end #hyperBess_j0prime


"""
	hyperBess_j0pprime(d, x)

Compute the 2nd derivative of the hyperspherical Bessel function of the first kind for
symmetrical states with respect to `x`. Symbolically, this is equal to:
	`- hyperBess_j(d, 0, x) / 2
		+ hyperBess_j(d, 2, x) / 2
		+ (d/2 - 1) * hyperBess_j(d, 1, x) / x`
This function is reserved for internal use for finding zeros of `hyperBess_j0prime`.
"""
#Float64 optimized implementation
function hyperBess_j0pprime(d::Float64, x::Float64)
	if d <= 0.0
		error("hyperBess_j0pprime: only positive values of d are supported")
	end

	baseorder::Float64 = muladd(0.5, d, -1.0)

	if x > abs(d) * 1e-12
		result = 0.5*x*(besselj(baseorder+2.0, x) - besselj(baseorder, x))
		result = muladd(baseorder, besselj(baseorder + 1.0, x), result)
		result *= x^(baseorder + 1.0)

	elseif 0.0 <= x #Do a first order approximation based on the definition
		result = muladd( 1.5 * x, x / (baseorder + 2.0), - 1.0) / d
		result /= 2.0^baseorder * gamma(baseorder + 1.0)

	else
		error("hyperBess_j0pprime: negative x outside of domain.")
	end

	return result
end#::Float64 #hyperBess_j0pprime


#Generic implementation
function hyperBess_j0pprime(d::Real, x::Tp) where Tp<:Number
	if d <= 0.0
		error("hyperBess_j0pprime: only positive values of d are supported")
	end
	baseorder = convert(Tp, d) / 2 - 1

	if x > abs(d) * 1e-12
		result = x*(besselj(real(baseorder)+2, x) - besselj(real(baseorder), x)) / 2
		result += baseorder * besselj(real(baseorder) + 1, x)
		result *= 2^baseorder * gamma(baseorder + 1) * x^(baseorder + 1)

	elseif 0 <= x #Do a first order approximation based on the definition
		result = (3 * x * x / (2*(baseorder + 2)) - 1) / d
		result /= 2^baseorder * gamma(baseorder + 1)

	else
		error("hyperBess_j0pprime: negative x outside of domain.")
	end

	return result
end #hyperBess_j0pprime


"""
	hyperBess_j0zeros(d, mmax)

Compute the first mmax zeros of the hyperspherical Bessel function of the first kind for
symmetrical states.
"""
#const MAX_NEWTON = convert(Int, )
#Float64 optimized version
function hyperBess_j0zeros( d::Float64, mmax::Int )
	if mmax <= 0
		error("hyperBess_j0zeros: mmax must be positive")
	end

	if d <= 0.0
		error("hyperBess_j0zeros: d must be positive")
	end

	baseorder::Float64 = muladd(0.5, d, -1.0)
	baseorderplus1::Float64 = baseorder + 1.0
	#empirical function that gets near the first zero
	guess::Float64 = sqrt((baseorder + 2.0)^2 - 1.0)

	result::Array{Float64,1} = Array{Float64}(mmax)
	converged::Bool = true

	#This function is not, exactly, hyperBess_j0, but it has the same zeros
	#and falls off more slowly at large x
	function f_fp(x::Float64)
		f::Float64 = besselj(baseorder, x)
		fp::Float64 = - besselj(baseorderplus1, x)

		return (f, fp)
	end

	for m in 1:mmax
		result[m], converged = findrootNewtonF64toF64(f_fp, guess,
				max_step_size=convert(Float64, pi/4), )
		if !converged
			error("hyperBess_j0zeros: failed to converge on zero $m with d=$d")
		end

		guess = result[m] + pi
	end


	return result
end#::Array{Float64,1} #hyperBess_j0zeros


#Generic implementation
function hyperBess_j0zeros( d::Real, mmax::Int, exampleval::Tp ) where Tp <: Real
	if mmax <= 0
		error("hyperBess_j0zeros: mmax must be positive")
	end

	if d <= 0
		error("hyperBess_j0zeros: d must be positive")
	end

	baseorder = convert(Tp, d)/2 - 1
	baseorderplus1 = convert(Tp, d)/2
	#empirical function that gets near the first zero
	guess = sqrt((baseorder + 2)^2 - 1)

	result = Array{Tp}(mmax)
	converged::Bool = true

	function f_fp(x::Tp)
		rootx = sqrt(x)
		f::Tp = besselj(baseorder, x) * rootx
		fp::Tp = - besselj(baseorderplus1, x) * rootx

		return (f, fp)
	end

	for m in 1:mmax
		result[m], converged = findrootNewton(f_fp, guess,
				max_step_size=convert(Tp, pi/4))
		if !converged
			error("hyperBess_j0zeros: failed to converge on zero $m with d=$d")
		end

		guess += pi
	end

	return result
end #hyperBess_j0zeros