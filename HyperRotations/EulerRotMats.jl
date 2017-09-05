"""
	matrixSONEuler( N, thetas )

	matrixSONEuler returns an N X N rotation matrix (a.k.a. an element of
	the group SO(N)), constructed from the angles stored in the 1-dimensional array,
	thetas. thetas must contain N * (N - 1) / 2 angles to uniquely specify an element
	of SO(N). The user is required to specify N because the dimension of the
	space is usually known, and the check for the validity of the length of
	thetas is faster this way.

	The angles in thetas are interpreted as defining a sequence of N - 1 unit
	vectors where every angle is 'polar.' For more information on how unit vectors
	are defined by polar angles see the documentation for the function
	UnitVecFromPolarAngles. The first angle fixes a 2-d unit vector, the
	next two a 3-d unit vector, etc all the way up to the N - 1 angles that
	fix an N-dimensional unit vector. It is this approach that makes the
	definition of the functions Jacobian_Euler and CanonicalizeEuler
	straightforward.

	Specifically, the rotation matrix, RotMat, is defined as the product
	of N - 1 rotation matrices:

	RotMat = R_2 * R_3 * ... * R_{N-1} * R_N

	R_i is a block-diagonal matrix, with an iXi rotation matrix in the
	upper left and the identity in the lower right. The rotation block
	of R_i is built with its last row as an i-dimensional unit vector,
	rhat_i, and its upper i - 1 rows are denoted thetahat_j. These rows
	are constructed as follows:

	thetahat_j = cos(theta_j) * rhat_j - sin(theta_j) * ehat_{j + 1}.

	As a concrete example, when N = 3, R_2 and R_3 are:

	`R_2 = [ cos(thetas[1])  -sin(thetas[1])  0;
		     sin(thetas[1])   cos(thetas[1])  0;
			 0                0               1 ]`
	`R_3 = [ cos(thetas[3])                   -sin(thetas[3])                    0             ;
	         cos(thetas[2]) * sin(thetas[3])   cos(thetas[2]) * cos(thetas[3])  -sin(thetas[2]);
	         sin(thetas[2]) * sin(thetas[3])   sin(thetas[2]) * cos(thetas[3])   cos(thetas[2]) ]`

	In the formalism used for the general rotate function in this package, the sequence
	of coordinate pairs that defines this function is:
	[ (1, 2),
	  (2, 3), (1, 2),
	  (3, 4), (2, 3), (1, 2),
	  ...
	  (N-1, N), (N-2, N-1), ..., (1, 2) ]
"""

function matrixSONEuler( N::Integer, thetas::Vector{Tp} ) where Tp<:Real
	lenAngs = length(thetas)

	if N * (N - 1) != 2 * lenAngs
		error("The number of dimensions (N) must match the number of angles (T):\n" *
			"\tN * (N - 1) = 2 * T")
	end

	if lenAngs < 1
		return Array{Tp}((0,0))
	end

	masterSinCosList = sincos(thetas)
	vecEnd = lenAngs
	RotMat = eye(Tp, N)
	R = Array{Tp}((N,N))
	rhat = Array{Tp}(N)
	thetahat = Array{Tp}(N)
	for row in N:-1:2
		fill!(R[row:N,1:N], zero(Tp))
		for i in (row+1):N
			R[i,i] = one(Tp)
		end

		vecStart = vecEnd - angsInVec + 1
		sincoses = masterSinCosList[vecStart:vecEnd]

		fill!(rhat[1:N], zero(Tp))
		rhat[1] = one(Tp)
		for i in 2:row
			s, c = sincoses[i - 1]

			thetahat = c .* rhat
			thetahat[i] = -s

			R[i-1,1:end] = thetahat

			rhat *= s
			rhat[i] = c
		end

		R[row,1:end] = rhat

		RotMat = R * RotMat

		vecEnd -= angsInVec
	end

	return RotMat
end


function anglesSONEuler!( RotMat::Matrix{Tp}, zeroTol::Tp = 8 * eps(Tp) ) where Tp<:Real
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
	vecEnd = iseven(N) ? div(N, 2) * (N - 1) : N * div((N - 1), 2)
	if vecEnd < N #The integer overflowed
		vecEnd = div(BigInt(N) * BigInt(N-1), 2)
	end

	thetas = Array{Tp}( vecEnd )
	R = Array{Tp}( size(RotMat) )
	rhat = Array{Tp}(N)
	thetahat = Array{Tp}(N)
	vecAngs = Array{Tp}(N-1)

	for row in N:-1:2
		angsInVec = row - 1
		vecStart = vecEnd - angsInVec + 1

		vecAngs[1:AngsInVec] = PolarAnglesFromUnitVec(RotMat[row, 1:row], zeroTol=zeroTol)
		thetas[vecStart:vecEnd] = vecAngs[1:AngsInVec]

		#construct the rotation matrix that projects the last row to the last unit vector
		# and its standard orthogonal compliment
		sincoses = sincos(vecAngs[1:AngsInVec])

		fill!(R[row:N,1:N], zero(Tp))
		for i in (row+1):N
			R[i,i] = one(Tp)
		end

		fill!(rhat[1:N], zero(Tp))
		rhat[1] = one(Tp)
		for i in 2:row
			s, c = sincoses[i - 1]

			thetahat = c .* rhat
			thetahat[i] = -s

			R[i-1,1:end] = thetahat

			rhat *= s
			rhat[i] = c
		end

		R[row,1:row] = RotMat[row,1:row]
		fill!(R[row,(row+1):N], zero(Tp))

		RowMat = RowMat * (R.')

		vecEnd -= angsInVec
	end

	return thetas
end


"""
	anglesSONEuler( RotMat, zeroTol=8 * eps() )

	anglesSONEuler is the inverse of matrixEuler. So, given an N X N rotation
	matrix, RotMat, this function returns the N * (N - 1) / 2 angles in the
	Euler parametrization of RotMat. The argument zeroTol sets the limit on what
	values are considered to be close enough to zero to treat as exactly zero in
	order to avoid the numerical instability intrinsic to the process of
	multiplying matrices by their numerical inverse. The anglesSONPolar! variant
	modifies RotMat in the process of the computation, turning it into an
	identity matrix, up to rounding errors.
	The anglesSONEuler! version uses the `RotMat` argument as scratch space.
"""

function anglesSONEuler( RotMat::Matrix{Tp}, zeroTol::Tp = 8 * eps(Tp) ) where Tp<:Real
	scratch = copy(RotMat)

	return anglesEuler!( scratch, zeroTol=zeroTol )
end