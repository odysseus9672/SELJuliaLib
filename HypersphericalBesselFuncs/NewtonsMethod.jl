
"""
    findrootNewtonF64toF64( func, funcder, guess, rtol=1e-12, abstol=1e-12, maxiter=100 )

Finds a root of `func` using Newton's method starting from `guess`, and signals
whether the algorithm converged to at least one of the tolerance goals within
`maxiter` loops. `funcder` must be the derivative of `func`. The iterative step is:

``x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}``

`rtol` is the relative tolerance goal for the search, and `abstol` is the absolute
tolerance goal, satisfying either of which is considered a successful search.
`func` and `funcder` must be type isomorphisms (maps from a subtype of Number to
the same subtype of Number), with that subtype fixed by the type of `guess`.
"""
#Float64 optimized version
function findrootNewtonF64toF64( func, funcder, guess::Float64,
	rtol::Float64=1e-12, abstol::Float64=1e-12, maxiter::Int=100 )
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

	converged::Bool = false

	oldx::Float64 = guess
	newx::Float64 = oldx - ((func(oldx) / funcder(oldx))::Float64)
	absdiff = abs(oldx - newx)

	iter::Int = 2
	while (absdiff > abstol || absdiff > rtol * abs(newx)) && iter <= maxiter
		oldx = newx
		newx = oldx - func(oldx) / funcder(oldx)
		if newx < 0.0
			newx = eps(Float64)
		end
		absdiff = abs(oldx - newx)

		iter += 1
	end #while (absdiff < abstol || absdiff < rtol * abs(newx)) && newxiter <= maxiter

	if iter <= maxiter
		converged = true
	end

	return (newx, converged)
end#::Tuple{Float64,Bool} #findrootNewtonF64toF64


#Generic version
function findrootNewton( func, funcder, guess::Tp,
						rtol::Real=1e-12, abstol::Real=1e-12, maxiter::Int=100 ) where Tp <: Number
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

	converged::Bool = false

	oldx = guess
	newx = oldx - ((func(oldx) / funcder(oldx))::Tp)
	absdiff = abs(oldx - newx)

	iter::Int = 2
	while (absdiff < abstol || absdiff < rtol * abs(newx)) && iter <= maxiter
		oldx = newx
		newx = oldx - func(oldx) / funcder(oldx)
		absdiff = abs(oldx - newx)

		iter += 1
	end #while (absdiff < abstol || absdiff < rtol * abs(newx)) && newxiter <= maxiter
	println( (absdiff < abstol, absdiff < rtol * abs(newx), iter <= maxiter) )
	if iter <= maxiter
		converged = true
	end

	return (newx, converged)
end #findzeroNewton