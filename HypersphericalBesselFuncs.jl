__precompile__()
module HypersphericalBesselFuncs

using SpecialFunctions

include("HypersphericalBesselFuncs/NewtonsMethod.jl")
include("HypersphericalBesselFuncs/hyperBess_j0.jl")
include("HypersphericalBesselFuncs/hyperBess_j.jl")
include("HypersphericalBesselFuncs/hyperBess_y0.jl")
include("HypersphericalBesselFuncs/hyperBess_y.jl")
include("HypersphericalBesselFuncs/hyperBess_i0.jl")
include("HypersphericalBesselFuncs/hyperBess_i.jl")
include("HypersphericalBesselFuncs/hyperBess_k0.jl")
include("HypersphericalBesselFuncs/hyperBess_k.jl")

export
	hyperBess_j0,
	hyperBess_j0G,
	hyperBess_j0prime,
	hyperBess_j0zeros,
	#hyperBess_j0primezeros,
	hyperBess_j,
	hyperBess_jprime,
	hyperBess_y0,
	hyperBess_y0prime,
	hyperBess_y0zeros,
	#hyperBess_y0primezeros,
	hyperBess_y,
	hyperBess_yprime,
	hyperBess_i0,
	hyperBess_i0prime,
	hyperBess_i,
	hyperBess_iprime,
	hyperBess_k0,
	hyperBess_k0prime,
	hyperBess_k,
	hyperBess_kprime


end #HypersphericalBesselFuncs