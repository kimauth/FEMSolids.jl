"""
    Neumann(f::Function, facetset::AbstractSet{FacetIndex}, nfacets::Int)

Represent a Neumann boundary condition.

# Arguments:
- `f`: Prescribed traction as a function `f(x,t)`, with spatial coordinate `x` and time `t`
- `facetset`: The facets on which the boundary condition acts
- `nfacets`: Number of facets per cell

Example:
```julia
grid = generate_grid(Quadrilateral, (5,5))

prescribed_traction(x,t) = Vec(zero(t), x[1]*t)
facetset = getfacetset(grid, "top")
nfacets = nfacets(getcells(grid, 1))

neumann_bc = Neumann(prescribed_traction, facetset, nfacets)
```
!!! note
    `Neumann` needs to be updated once for each time step with [`update!`](@ref) and
    once for each cell with [`update_cell!`](@ref). 
"""
mutable struct Neumann{F, FS <: AbstractSet{FacetIndex}}
    const f::F
    const facetset::FS
    const nfacets::Int
    time::Float64
    current_cellid::Int
end

function Neumann(f::Function, facetset::AbstractSet{FacetIndex}, nfacets)
    return Neumann(f, facetset, nfacets, 0.0, 0)
end

"""
    update!(bc::Neumann, time::Float64=0.0)

Update the Neumann boundary condition `bc` for the current `time`.
"""
function Ferrite.update!(bc::Neumann, time::Float64 = 0.0)
    bc.time = time
    return bc
end

"""
    update_cell!(bc::Neumann, cellid::Int)

Update the Neumann boundary condition `bc` with the current `cellid`.
"""
function update_cell!(bc::Neumann, cellid::Int)
    bc.current_cellid = cellid
    return bc
end
