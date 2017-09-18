
"""
    findrootNewtonF64toF64( f_fprime, guess, rtol=1e-12, abstol=1e-12, maxiter=100 )

Finds a root of `f` using Newton's method starting from `guess`, and signals
whether the algorithm converged to at least one of the tolerance goals within
`maxiter` loops. `f_fprime` must return both the function, `f`, and its derivative,
`fprime`. The iterative step is:

``x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}``

`rtol` is the relative tolerance goal for the search, and `abstol` is the absolute
tolerance goal, satisfying either of which is considered a successful search.
`func` and `funcder` must be type isomorphisms (maps from a subtype of Number to
the same subtype of Number), with that subtype fixed by the type of `guess`.
This implementation has been tuned to be optimized for finding zeros of Bessel
functions.
"""
#Float64 optimized version
function findrootNewtonF64toF64( f_fprime, guess::Float64;
	rtol::Float64=16*eps(Float64), abstol::Float64=eps(Float64)/16, maxiter::Int=100,
	max_step_size::Float64=0.25 * pi )
	#make sure tolerances are positive
	if rtol <= 0.0
		error("findrootNewton: rtol must be a positive number")
	end

	if abstol <= 0.0
		error("findrootNewton: abstol must be a positive number")
	end

	if maxiter <= 0
		error("findrootNewton: maxiter must be a positive integer")
	end

	if max_step_size <= 0.0
		error("findrootNewton: max_step_size must be a positive number")
	end

	converged::Bool = false

	oldx::Float64 = guess
	f::Float64, fp::Float64 = f_fprime(oldx)
	Newton_step::Float64 = f/fp
	absdiff::Float64 = min(abs(Newton_step), max_step_size)
	newx::Float64 = oldx - sign(Newton_step) * absdiff


	iter::Int = 2
	while absdiff > abstol && absdiff > rtol * abs(newx) && iter <= maxiter
		oldx = newx
		f, fp = f_fprime(oldx)
		Newton_step = f/fp
		absdiff = min(abs(Newton_step), max_step_size)
		newx = oldx - sign(Newton_step) * absdiff

		if newx < 0.0
			newx = eps(Float64)
			absdiff = abs(newx - oldx)
		end

		iter += 1
	end #while (absdiff < abstol || absdiff < rtol * abs(newx)) && newxiter <= maxiter

	if iter <= maxiter
		converged = true
	else
		println( (guess, oldx, newx) )
	end

	return (newx, converged)
end#::Tuple{Float64,Bool} #findrootNewtonF64toF64


#Generic version
function findrootNewton( f_fprime, guess::Tp;
						rtol::Tp2=16 * eps(Tp2), abstol::Tp3=eps(Tp3)/16,
						maxiter::Int=100,
						max_step_size = convert(Tp, pi/4) ) where {Tp <: Number, Tp2 <: Real, Tp3 <: Real}
	#make sure tolerances are positive
	if rtol <= 0
		error("findrootNewton: rtol must be a positive number")
	end

	if abstol <= 0
		error("findrootNewton: abstol must be a positive number")
	end

	if maxiter <= 0
		error("findrootNewton: maxiter must be a positive integer")
	end

	if max_step_size <= 0
		error("findrootNewton: max_step_size must be a positive number")
	end

	converged::Bool = false

	oldx = guess
	f, fp = f_fprime(oldx)
	Newton_step = f/fp
	absdiff = min(abs(Newton_step), max_step_size)
	newx = oldx - sign(Newton_step) * absdiff


	iter::Int = 2
	while absdiff > abstol && absdiff > rtol * abs(newx) && iter <= maxiter
		oldx = newx
		f, fp = f_fprime(oldx)
		Newton_step = f/fp
		absdiff = min(abs(Newton_step), max_step_size)
		newx = oldx - sign(Newton_step) * absdiff
		absdiff = abs(oldx - newx)

		iter += 1
	end #while (absdiff < abstol || absdiff < rtol * abs(newx)) && newxiter <= maxiter
	println( (absdiff < abstol, absdiff < rtol * abs(newx), iter <= maxiter) )
	if iter <= maxiter
		converged = true
	end

	return (newx, converged)
end #findzeroNewton