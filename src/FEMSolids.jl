module FEMSolids
using Ferrite, Tensors

include("linear_elasticity.jl")
include("Neumann.jl")
include("element_routine.jl")
include("adaptivity.jl")

export LinearElasticity
export Neumann, update_cell!
export bulk_modulus, shear_modulus
export element_routine!, Primal
export refine, transfer_solution

end
