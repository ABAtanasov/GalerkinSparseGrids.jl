# Terms Used

"the interval"				The unit interval [0,1]

"interval"					A nonzero sub-interval of [0,1], possibly all of it
("sub-interval")			Always dyadic in our case, of the form (i/2^n, (i+1)/2^n)

"interval product"			A cartesian product of D sub-intervals in D dimensions

"resolution"				An interval of the form (i/2^n, (i+1)/2^n) has resolution n
							When n = 0 this involves functions that can be nonzero on all of [0,1]
							In order to count hat basis functions that do not vanish on the boundaries 
							we label them with resolution(s) '-1' (and '-2' for hat basis) for simplicity
							In the Galerkin basis, resolution '-1' corresponds to Legendre polynomials on [0,1]
							

"level"						The level of resolution. This is just n above, when n >= 0, and 
							
"cell"						For a given level l of resolution, a specifier for one of the above sub-intervals

"mode"						In the Galerkin Scheme of degree k, 
							each (level, cell) specifies a sub-interval on which we have k basis functions 
							(polynomials up to degree p = k-1). 
							mode uniquely determines one of the k such basis functions


"multi-level"				A vector of D numbers specifying the resolution along each of D dimensions
variable name: "levels"

"multi-cell"				A vector of D numbers specifying a specific interval product given a multilevel
variable name: "cells"

"multi-mode"				A vector of D numbers specifying a specific product of galerkin basis functions for 
variable name: "modes"		each interval associated with the interval product obtained from (levels, cells)
							

"block"						The space of basis functions associated to a multilevel