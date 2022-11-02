"""
    Neumann(f::Function, faceset::Set{FaceIndex}, nfaces::Int)

Represent a Neumann boundary condition.

# Arguments:
- `f`: Prescribed traction as a function `f(x,t)`, with spatial coordinate `x` and time `t`
- `faceset`: The faces on which the boundary condition acts
- `nfaces`: Number of faces per cell

Example:
```julia
grid = generate_grid(Quadrilateral, (5,5))

prescribed_traction(x,t) = Vec(zero(t), x[1]*t)
faceset = getfaceset(grid, "top")
nfaces = nfaces(getcells(grid, 1))

neumann_bc = Neumann(prescribed_traction, faceset, nfaces)
```
!!! note
    `Neumann` needs to be updated once for each time step with [`update!`](@ref) and
    once for each cell with [`update_cell!`](@ref). 
"""
struct Neumann{F}
    f::F
    faceset::Set{FaceIndex}
    nfaces::Int
    time::Ferrite.ScalarWrapper{Float64}
    current_cellid::Ferrite.ScalarWrapper{Int}
end

function Neumann(f::Function, faceset::Set{FaceIndex}, nfaces)
    return Neumann(f, faceset, nfaces, Ferrite.ScalarWrapper(0.0), Ferrite.ScalarWrapper(0))
end

"""
    update!(bc::Neumann, time::Float64=0.0)

Update the Neumann boundary condition `bc` for the current `time`.
"""
function Ferrite.update!(bc::Neumann, time::Float64 = 0.0)
    bc.time[] = time
    return bc
end

"""
    update_cell!(bc::Neumann, cellid::Int)

Update the Neumann boundary condition `bc` with the current `cellid`.
"""
function update_cell!(bc::Neumann, cellid::Int)
    bc.current_cellid[] = cellid
    return bc
end
