module FEMSolids

using Reexport
@reexport using Ferrite

include("linear_elasticity.jl")
include("Neumann.jl")
include("element_routine.jl")
include("adaptivity.jl")

export LinearElasticity
export Neumann, update_cell!
export bulk_modulus, shear_modulus
export element_routine!, Primal
export refine, linear_to_quadratic

end
