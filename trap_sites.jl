
include("mcclean_isotherm_conc_dist.jl")
using McCleanIsotherm

""" This file describes the trap site occupation and the
redistribution of the trap site probabilities

The spline function which describes the concentration dependence comes
from McCleanIsotherm. This function will be used to define the solute
concentration in each dislocation segment along the line.
 """
