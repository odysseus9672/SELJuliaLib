"""
	hypotVec( invec )

	Computes the Euclidean norm of its vector argument using a
	generalization of the algorithm used to compute the hypotenuse of
	2-dimensional vectors in the built in function hypot(x, y), and
	adjusted so that it should also handle lists which obey a triangle inequality
	(ie no single element is unusually large). It's recursive definition is:

	`hypot(vec) = ( hypot(big half of vec) *
					sqrt( 1 + sum((small half of vec)^2) / sum((big half of vec)^2) ) )
	hypot(vec of length 1) = abs(vec[1])`

	(small half of vec) is the first half of vec, rounded down if its length is odd,
	after vec's elements have been sorted in increasing absolute value.
	(big half of vec) is the elements from vec not in (small half of vec)."""

function hypotVec( invec::Vector{Tp} ) where Tp<:Real
	if any(isnan(invec))
		return NaN
	end

	if any(isinf(invec))
		return Inf
	end

	veclen::Int = length(invec)
	result::Tp
	if veclen > 2
		sortedvec = sort(invec, by=abs)

		#Calculate the size of the scratch area needed
		lenScratch::Int = 1 + ceil(log2(veclen))
		Scratch = zeros(Tp, lenScratch)

		#Scale the sorted vec to control possible overflow
		biggestInput::Tp = abs(sortedvec[veclen])
		sortedvec /= biggestInput

		#Produce the intermediate sums
		stoppoint::Int = 0
		iScratch::Int = 1
		while iScratch <= lenScratch
			startpoint = stoppoint + 1
			compstodigest::Int = veclen - startpoint + 1
			stoppoint = max(compstodigest >> 1, 1)
			stoppoint += startpoint - 1

			Scratch[iScratch] = sum_kbn(sortedvec[startpoint:stoppoint].^2)

			iScratch += 1
		end

		#Digest the scratch into results
		result = 1
		iScratch = 1
		while iScratch < lenScratch
			bigpart = sum_kbn(Scratch[iScratch+1:lenScratch])

			if Scratch[iScratch] > 0
				result *= 1 + Scratch[iScratch] / bigpart
			end

			iScratch += 1
		end

		result = biggestInput * sqrt(result)

	elseif veclen == 2
		result = hypot(invec[1], invec[2])
	elseif veclen == 1
		result = invec[1]
	else
		result = NaN
	end

	return result
end
