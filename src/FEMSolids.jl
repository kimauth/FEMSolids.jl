module FEMSolids

using Reexport
@reexport using Ferrite

include("linear_elasticity.jl")
include("element_routine.jl")

export LinearElasticity
export element_routine!, Primal

end
