"""
	matrixSONTaitBryan( N, thetas )

	Returns an SO(N) rotation matrix that is a generalization of the Tait-Bryan
	system in SO(3).

	In the formalism used for the general rotate function in this package, the sequence
	of coordinate pairs that defines this function is:
	[ (1, 2), (1, 3), ..., (1, N),
	  (2, 3), (2, 4), ..., (2, N),
	  ...
	  (N-2, N-1), (N-2, N),
	  (N-1, N) ]
"""

function matrixSONTaitBryan( N::Integer, thetas::Vector{Tp} ) where Tp<:Real
	lenAngs = length(thetas)

	if N * (N - 1) != 2 * lenAngs
		error("The number of dimensions (N) must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	if lenAngs < 1
		return Array{Tp}((0,0))
	end

	coordPairs = [ (1, j) for j in 2:N ]
	for i in 2:(N-1)
		coordPairs = vcat(coordPairs, [ (i, j) for j in (i+1):N ])
	end

	RotMat = eye(Tp, N)
	rotateSeq!(RotMat, 2, thetas, coordPairs)

	return RotMat
end


#=
export anglesSONTaitBryan!

function anglesSONTaitBryan!{Tp<:FloatingPoint}( RotMat::Matrix{Tp}, zeroTol::Tp = 8 * eps(Tp) )
	#Test if the matrix is square
	N, M = size(RotMat)
	if N != M
		DomainError()
	end

	#Test if the matrix is orthogonal, to within NumerTol
	if any(abs( RotMat.' * RotMat - convert(Array{Tp}, eye(size(RotMat))) )) < zeroTol
		DomainError()
	end

	#Minimize integer overflow in calculating N * (N - 1) / 2 by testing for N even
	numthetas = iseven(N) ? div(N, 2) * (N - 1) : N * div((N - 1), 2)
	if numthetas < N #The integer overflowed
		numthetas = div(BigInt(N) * BigInt(N-1), 2)
	end

	thetas = Array{Tp}( numthetas )

end

export anglesSONTaitBryan{Tp<:FloatingPoint}( RotMat::Matrix{Tp}, zeroTol::Tp = 8 * eps(Tp) )
	scratch = copy(RotMat)

	return anglesSONTaitBryan!( scratch, zeroTol=zeroTol )
end
=#