"""
	rotateSO3ypr!(v, yaw, pitch, roll)

	rotateSO3ypr! rotates the vector stored in its first argument
	according to the Tait-Bryan angles using to the ZYX convention.
"""

function rotateSO3ypr!(v::Vector{Tp}, yaw::Tp, pitch::Tp, roll::Tp) where Tp<:Real
	if length(v) != 3 || length(axis) != 3
		DomainError()
	end

	sy, cy = sincos(yaw)
	sp, cp = sincos(pitch)
	sr, cr = sincos(roll)

	srmy = sin(roll - yaw)
	crmy = cos(roll - yaw)
	cvp = coversin(pitch)

	x, y, z = v[1], v[2], v[3]
	v[1] = ( cp * cy) * x + (-sr * cy * cvp + srmy) * y + (-cr * cy * cvp + crmy) * z
	v[2] = ( cp * sy) * x + (-sr * sy * cvp + crmy) * y + (-cr * sy * cvp - srmy) * z
	v[3] = ( sp) * x + ( cp * sr) * y + ( cp * cr) * z

	return
end


"""
	rotateSO3ypr(v, yaw, pitch, roll)

	rotateSO3ypr returns a copy of the vector stored in its first argument
	rotated according to the Euler angles using to the ZYX convention.
"""

function rotateSO3ypr(v::Vector{Tp}, yaw::Tp, pitch::Tp, roll::Tp) where Tp<:Real
	res = copy(v)
	rotateSO3ypr!(res, yaw, pitch, roll)
	return res
end


"""
	rotateSO3Euler!(v, alpha, beta, gamma)

	rotateSO3Euler! rotates the vector stored in its first argument
	according to the Euler angles using to the ZXZ convention.
"""

function rotateSO3Euler!(v::Vector{Tp}, alpha::Tp, beta::Tp, gamma::Tp) where Tp<:Real
	if length(v) != 3 || length(axis) != 3
		DomainError()
	end

	sa, ca = sincos(alpha)
	sb, cb = sincos(beta)
	sg, cg = sincos(gamma)

	capg = cos(alpha + gamma)
	sapg = sin(alpha + gamma)
	vb = coversin(beta)

	x, y, z = v[1], v[2], v[3]
	v[1] = ( sa * sg * vb + capg) * x + ( sa * cg * vb - sapg) * y + ( sa * sb) * z
	v[2] = (-ca * sg * vb + sapg) * x + (-ca * cg * vb - capg) * y + (-ca * sb) * z
	v[3] = ( sb * sg) * x + ( sb * cg) * y + ( cb) * z

	return
end


"""
	rotateSO3Euler(v, alpha, beta, gamma)

	rotateSO3Euler returns a copy of the vector stored in its first argument
	rotated according to the Euler angles using to the ZXZ convention.
"""

function rotateSO3Euler(v::Vector{Tp}, alpha::Tp, beta::Tp, gamma::Tp) where Tp<:Real
	res = copy(v)
	rotateSO3Euler!(res, alpha, beta, gamma)
	return res
end


"""
	rotateSO3AxisAng!(v, axis, theta)

	rotateSO3AxisAng! rotates the vector stored in its first argument
	counter clock-wise around the axis in its second argument (does not need to be
	normalized) by the angle theta.
"""

function rotateSO3AxisAng!(v::Vector{Tp}, axis::Vector{Tp}, theta::Tp) where Tp<:Real
	if length(v) != 3 || length(axis) != 3
		DomainError()
	end

	st, ct = sincos(theta)
	versint = versin(theta)
	uvec = NormalizeVec(axis)

	v[1:3] = ct * v + st * cross(axis, v) + (versint * dot(v, axis)) * v
end


"""
	rotateSO3AxisAng(v, axis, theta)

	rotateSO3AxisAng returns a copy of the vector stored in its first argument
	rotated counter clock-wise around the axis in its second argument (does not need to be
	normalized) by the angle theta.
"""

function rotateSO3AxisAng(v::Vector{Tp}, axis::Vector{Tp}, theta::Tp) where Tp<:Real
	result = copy(v)
	rotate3VecAxisAng!(result, axis, theta)
	return result
end