# SELJuliaLibs

This Github contains some utility libraries I made for entertainment. Though
they are intended to be useful for the purposes described in their
internal documentation, I make no warranty that they'll be useful for any
particular purpose.

The HyperRotations module is intended to contain a number of functions
useful for producing N by N rotation matrices, and rotating vectors in
N-dimensions.

The HypersphericalBesselFuncs module contains Bessel functions generalized
for `d`-dimensions. Primarily, this is a wrapper around functions in the
Julia `SpecialFunctions` package, and not all special cases have been debugged,
and therefore might not be handled correctly, but it should work.

SEL
2017-09-04