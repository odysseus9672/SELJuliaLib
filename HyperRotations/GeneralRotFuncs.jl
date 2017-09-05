using ArrayViews

function rotateSeq!(v::Vector{Tp}, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	if length(thetas) != size(coordPairs)[1]
		DimensionMismatch()
	end

	for i in length(thetas):-1:1
		c, s = sincos(thetas[i])
		fromidx = coordPairs[i,1]
		toidx = coordPairs[i,2]
		steps = toidx - fromidx
		oddness = isodd(steps)

		if steps == 0
			error("Rotations can only be defined as mixing distinct components")
		end

		if (oddness && steps < 0) || (!oddness && steps > 0)
			s *= -one(Tp)
		end

		v0, v1 = v[fromidx], v[toidx]

		v[fromidx], v[toidx] = v0 * c - v1 * s, v1 * c + v0 * s
	end
end

function rotateSeq!(v::Matrix{Tp}, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	if length(thetas) != size(coordPairs)[1]
		DimensionMismatch()
	end

	numvecs = size(v)[2]
	for i in length(thetas):-1:1
		c, s = sincos(thetas[i])
		fromidx = coordPairs[i,1]
		toidx = coordPairs[i,2]
		steps = toidx - fromidx
		oddness = isodd(steps)

		if steps == 0
			error("Rotations can only be defined as mixing distinct components")
		end

		if (oddness && steps < 0) || (!oddness && steps > 0)
			s *= -one(Tp)
		end

		v0, v1 = v[fromidx,1:end], v[toidx,1:end]

		v[fromidx,1:end], v[toidx,1:end] = v0 * c - v1 * s, v1 * c + v0 * s
	end
end

function rotateSeq!(v::Array{Tp}, vecidx::Int, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	if length(thetas) != size(coordPairs)[1]
		DimensionMismatch()
	end

	r = length(size(v)) #ie rank of v as a tensor
	if vecidx > length(size(v))
		error("vecidx must be a valid index to size(v)")
	end

	pfx = [ Colon() for i in 1:(vecidx - 1) ]
	sfx = [ Colon() for i in (vecidx + 1):r ]

	for i in length(thetas):-1:1
		c, s = sincos(thetas[i])
		fromidx = coordPairs[i,1]
		toidx = coordPairs[i,2]
		steps = toidx - fromidx
		oddness = isodd(steps)

		if steps == 0
			error("Rotations can only be defined as mixing distinct components")
		end

		if (oddness && steps < 0) || (!oddness && steps > 0)
			s *= -one(Tp)
		end

		v0 = view(v, tuple(vcat(pfx, fromidx, sfx)...))
		v1 = view(v, tuple(vcat(pfx, toidx, sfx)...))

		v0, v1 = v0 * c - v1 * s,  v0 * s + v1 * c
	end
end


"""
	rotateSeq(v::Vector{Tp}, thetas::Vector{Tp}, coordPairs::Array{Tp, 2})
	rotateSeq(v::Matrix{Tp}, thetas::Vector{Tp}, coordPairs::Array{Tp, 2})
	rotateSeq(v::Array{Tp}, vecidx::Int, thetas::Vector{Tp}, coordPairs::Array{Tp, 2})

	Applies a sequence of elemental rotations to the vector(s) stored in
	its first argument. If the first argument is a matrix, then the vectors
	are assumed to be stored in the rows. Optionally, an integer second
	argument can be added to specify along which axis the vectors are stored.
	The rotations are applied in reverse order, as though the elemental
	rotation matrices were constructed and multiplied in the same order as
	the angles stored in thetas and then right multiplied on the vector, ie

	vec_rotated = RotMatrix * vec.

	Any general rotation can be built up using multiple different sequences
	of elemental rotations.

	An elemental rotation is defined by an angle and two
	integers specifying the components it mixes. Note that this is unlike 3-d
	where specifying the axis not-rotated is sufficient, becuase the invariant
	space here has dimension N - 2. For an N-dimensional elemental rotation the
	detailed definition is as follows:
	Let t be the angle of rotation, and let i and j be the components being
	mixed with i designated the 'from' component, and j the 'to.' If j = i + 1
	(ie the to component is immediately after the from component) then the rotation
	has the standard 2-d rotation matrix form:
	v[i], v[j] = cos(t) * v[i] - sin(t) * v[j], sin(t) * v[i] + cos(t) * v[j].
	Otherwise, the sign of t is flipped the number of times adjacent components
	would need to be flipped in order arrange from the from component to
	immediately precede the to component. This is why, for example, the 3-d
	rotation about the y-axis (a rotation from x to z) has the form:
	[  cos(t)  0  sin(t);
	   0       1  0     ;
	  -sin(t)  0  cos(t) ]
	where the minus sign is on the opposite side of the diagonal compared to
	rotations about the x and z axes."""

function rotateSeq(v::Vector{Tp}, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	result = copy(v)

	rotate!(result, thetas, coordPairs)

	return result
end

function rotateSeq(v::Matrix{Tp}, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	result = copy(v)

	rotate!(result, thetas, coordPairs)

	return result
end

function rotateSeq(v::Array{Tp}, vecidx::Int, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	result = copy(v)

	rotate!(result, vecidx, thetas, coordPairs)

	return result
end


"""
	matrixRotateSeq(N, thetas, coordPairs)

	matrixRotateSeq constructs the rotation matrix built up by a sequence of
	elemental rotations. See the documentation of the function `rotateSeq` for
	details.
"""

function matrixRotateSeq(N::Integer, thetas::Vector{Tp}, coordPairs::Array{Tp, 2}) where Tp<:Real
	if length(thetas) != size(coordPairs)[1]
		DimensionMismatch()
	end

	Rot = eye(Tp, N, N)

	rotateSeq!(Rot, thetas, coordPairs)

	return Rot
end
