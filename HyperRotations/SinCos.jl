#Julia translation of SinCos.go from: http://golang.org/src/pkg/math/sincos.go

const PI4A = 0x1.921fb4p-1
const PI4B = 0x1.4442dp-25
const PI4C = 0x1.8469898cc5170p-49
const M4PI = 0x1.45f306dc9c882p+0
const fallbackangle = 0x1p20 #Angle for which fallback to built in functions occurs

#sin coefficients
const sincoeffs = [ 0x1.5d8fd1fd19ccdp-33,
				   -0x1.ae5e5a9291f5dp-26,
				    0x1.71de3567d48a1p-19,
				   -0x1.a01a019bfdf03p-13,
				    0x1.111111110f7d0p-7,
				   -0x1.5555555555548p-3,
				    0x1.0p0 ]

const coscoeffs = [-0x1.8fa49a0861a9bp-37,
				    0x1.1ee9d7b4e3f05p-29,
				   -0x1.27e4f7eac4bc6p-22,
				    0x1.a01a019c844f5p-16,
				   -0x1.6c16c16c14f91p-10,
				    0x1.555555555554bp-5,
				   -0x1.0p-1 ]


"""
	sincos( theta )

	Calculates both the sine and cosine of the angle `theta`, in radians.
	If `theta` is a floating point number that can be represented without
	loss of accuracy as a `Float64`, then the calculation should be more
	efficient than calling sin and cos separately. Otherwise, sin and cos are
	called separately. The two values are returned in a tuple, sine first and
	cosine second.
"""

function sincos( theta::Tp ) where Tp<:Real
	return (sin(theta), cos(theta))
end

#Optimized versions
function sincos( theta::Float16 )
	return convert(Float16, _sincos64(convert(Float64, theta)))
end
function sincos( theta::Float32 );
	return convert(Float32, _sincos64(convert(Float64, theta)))
end
function sincos( theta::Float64 ); _sincos64(theta); end

function _sincos64( theta::Float64 )
	signofsin::Int32 = 0
	signofcos::Int32 = 0
	y::Float64 = 0.0
	z::Float64 = 0.0
	zz::Float64 = 0.0
	sinval::Float64
	cosval::Float64
	j::Int64

	#Handle special cases
	if theta == 0.0
		return zero(Tp), one(Tp)
	elseif isnan(theta) || isinf(theta)
		return convert(Tp, NaN), convert(Tp, NaN)
	elseif abs(theta) > fallbackangle
		return sin(theta), cos(theta)
	end

	signofsin = 0
	if theta > 0.0
		signofsin = 0
	else
		signofsin = 1
		theta = -theta
	end

	#integer part of x/(Pi/4), as integer for tests on the phase angle
	j = convert(Int64, trunc(M4PI * theta))

	# map zeros to origin
	if j&1 == 1
		j += 1
	end

	# integer part of x/(Pi/4), as float
	y = convert(Float64, j)

	#octant modulo 2Pi radians (360 degrees)
	j &= 7

	#reflect in x axis
	if j > 3
		j -= 4
		signofsin $= 1
		signofcos $= 1
	end
	if j > 1
		signofcos $= 1
	end

	#Extended precision modular arithmetic
	z = ((theta - y*PI4A) - y*PI4B) - y*PI4C
	zz = z * z

	cosval = coscoeffs[1]
	sinval = sincoeffs[1]

	for i = 2:7
		cosval = cosval * zz + coscoeffs[i]
		sinval = sinval * zz + sincoeffs[i]
	end

	cosval = cosval * zz + 1.0
	sinval = sinval * z + 0.0 #maintains FMA form

	if j == 1 || j == 2
		sinval, cosval = cosval, sinval
	end

	if signofcos == 1
		cosval = -cosval
	end
	if signofsin == 1
		sinval = -sinval
	end

	return (sinval, cosval)
end


"""
	sinver( theta )

	Calculates both the sine and versine of the angle `theta`, in radians.
	If `theta` is a floating point number that can be represented without loss of
	accuracy as a `Float64`, then the calculation should be more efficient than
	calling sin and versine separately. Otherwise, sin and versin are called
	separately. The two values are returned in a tuple, sine first and versine
	second.
"""

function sinver( theta::Tp ) where Tp<:Real
	return (sin(theta), versin(theta))
end

#Optimized versions
function sinver( theta::Float16 )
	return convert(Float16, _sinver64(convert(Float64, theta)))
end
function sinver( theta::Float32 );
	return convert(Float32, _sinver64(convert(Float64, theta)))
end
function sinver( theta::Float64 ); _sinver(theta); end

function _sinver64( theta::Float64 )
	signofsin::Int = 0
	signofcos::Int = 0
	y::Float64 = 0.0
	z::Float64 = 0.0
	zz::Float64 = 0.0
	sinval::Float64
	verval::Float64
	j::Int

	#Handle special cases
	if theta == 0.0
		return zero(Tp), one(Tp)
	elseif isnan(theta) || isinf(theta)
		return convert(Tp, NaN), convert(Tp, NaN)
	elseif abs(theta) > fallbackangle
		return sin(theta), verval(theta)
	end

	signofsin = 0
	if theta > 0.0
		signofsin = 0
	else
		signofsin = 1
		theta = -theta
	end

	#integer part of x/(Pi/4), as integer for tests on the phase angle
	j = convert(Int64, trunc(M4PI * theta))

	# map zeros to origin
	if j&1 == 1
		j += 1
	end

	# integer part of x/(Pi/4), as float
	y = convert(Float64, j)

	#octant modulo 2Pi radians (360 degrees)
	j &= 7

	#reflect in x axis
	if j > 3
		j -= 4
		signofsin $= 1
		signofcos $= 1
	end
	if j > 1
		signofcos $= 1
	end

	#Extended precision modular arithmetic
	z = ((theta - y*PI4A) - y*PI4B) - y*PI4C
	zz = z * z

	verval = coscoeffs[1]
	sinval = sincoeffs[1]

	for i = 2:7
		verval = verval * zz + coscoeffs[i]
		sinval = sinval * zz + sincoeffs[i]
	end

	verval *= -zz
	sinval *= z

	if j == 1 || j == 2
		sinval, verval = 1.0 - verval, 1.0 - sinval
	end

	if signofcos == 1
		verval = 2.0 - verval
	end
	if signofsin == 1
		sinval = -sinval
	end

	return (sinval, verval)
end
