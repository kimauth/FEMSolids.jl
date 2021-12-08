using FEMSolids
using Test

@testset "FEMSolids.jl" begin
    include("test_element_routine.jl")
    include("test_adaptivity.jl")
end
