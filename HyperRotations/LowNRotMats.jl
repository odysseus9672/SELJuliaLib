"""
	matrixSO2(theta)

	matrixSO2 produces the 2-d rotation matrix. If this function
	returns ``R``, and if ``v_{\\mathrm{new}} = R v_{\\mathrm{old}}``,
	then ``v_{\\mathrm{new}} is vold rotated counter-clockwise by theta radians.
"""

function matrixSO2(theta::Tp) where Tp<:Real
	s, c = sincos(theta)
	return Array{Tp}[ c -s; s c ]
end


"""
	matrixSO3ypr(yaw, pitch, roll)
	matrixSO3ypr returns the 3-d rotation matrix defined by the standard
	yaw-pitch-roll description.
"""

function matrixSO3ypr(yaw::Tp, pitch::Tp, roll::Tp) where Tp<:Real
	sy, cy = sincos(yaw)
	sp, cp = sincos(pitch)
	sr, cr = sincos(roll)

	srmy = sin(roll - yaw)
	crmy = cos(roll - yaw)
	cvp = coversin(pitch)

	return Tp[ ( cp * cy )  (-sr * cy * cvp + srmy )  (-cr * cy * cvp + crmy ) ;
			   ( cp * sy )  (-sr * sy * cvp + crmy )  (-cr * cy * cvp - srmy ) ;
			   (-sp )       ( cp * sr )               ( cp * cr )                ]

	#Naive version:
	#=return Tp[ (cp * cy)  (sp * sr * cy - cr * sy)  (sr * sy + sp * cr * cy) ;
			   (cp * sy)  (sp * sr * sy + cr * cy)  (sp * cr * sy - sr * cy) ;
			   (-sp)      (cp * sr)                 (cp * cr)                 ]=#
end


"""
	matrixSO3Euler(alpha, beta, gamma)

	matrixSO3Euler returns the 3d rotation matrix defined by the
	Euler angles alpha, beta, and gamma. If Rx(t) is a rotation around the
	x-axis by angle t, then this function returns:
	Ry(gamma) * Rx(beta) * Ry(alpha).
	Other Euler angle rotations are related to this one by various permutations
	of the axes.
"""

function matrixSO3Euler(alpha::Tp, beta::Tp, gamma::Tp) where Tp<:Real
	sa, ca = sincos(alpha)
	sb, cb = sincos(beta)
	sg, cg = sincos(gamma)

	capg = cos(alpha + gamma)
	sapg = sin(alpha + gamma)
	vb = versin(beta)
	#Not sure under what conditions this is more accurate than the naive version
	return Tp[ ( sa * sg * vb + capg )  ( sa * cg * vb - sapg )  ( sa * sb ) ;
			   (-ca * sg * vb + sapg )  (-ca * cg * vb + capg )  (-ca * sb ) ;
			   ( sb * sg )              ( sb * cg )              ( cb )       ]

	#Naive version:
	#=return Tp[ (ca * cg - sa * cb * sg)  (-ca * sg - sa * cb * cg)  (sa * sb)  ;
			   (ca * cb * sg + sa * cg)  (ca * cb * cg - sa * sg)   (-ca * sb) ;
			   (sb * sg)                 (sb * cg)                  (cb)        ]=#
end


"""
	matrixSO3AxisAng(axis, theta)

	matrixSO3AxisAng produces the matrix that represents a counter-clockwise
	rotation about the axis specified by the vector axis (doesn't need to be a
	unit vector) by the angle theta.
"""

function matrixSO3AxisAng(axis::Vector{Tp}, theta::Tp) where Tp<:Real
	if length(axis) != 3
		DomainError()
	end

	st, ct = sincos(theta)
	versint = versin(theta)
	uvec = NormalizeVec(axis)
	ux = uvec[1]
	uy = uvec[2]
	uz = uvec[3]

	vx = ux * versint
	vy = uy * versint
	vz = uz * versint

	sx = ux * st
	sy = uy * st
	sz = uz * st

	return [ ( ct + ux*vx)  (-sz + uy*vx)  ( sy + uz*vx);
			 ( sz + ux*vy)  ( ct + uy*vy)  (-sx + uz*vy);
			 (-sy + ux*vz)  ( sx + uy*vz)  ( ct + uz*vz) ]
end
