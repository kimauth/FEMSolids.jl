"""
    LinearElasticity{dim}(;G, K) where dim

Construct a `LinearElasticity` material based on the shear modulus `G` and the bulk modulus `K`.
"""
struct LinearElasticity{T,dim,M}
    G::T
    K::T
    Eᵉ::SymmetricTensor{4,dim,T,M}
end

function LinearElasticity{dim}(;G::Float64, K::Float64) where dim
    λ = K - 2/3*G
    δ(i,j) = i==j ? 1.0 : 0.0
    f(i,j,k,l) = λ*δ(i,j)*δ(k,l) + G*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    E = SymmetricTensor{4, dim}(f)
    return LinearElasticity(G, K, E)
end

bulk_modulus(;E, ν) = E / 3.0(1.0-2.0ν)
shear_modulus(;E, ν) = E / 2.0(1.0+ν)
