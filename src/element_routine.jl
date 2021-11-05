
"""
    Primal

This struct is used for dispatching which `element_routine!` is used when you add your own element routines.
`Primal` refers to the primal format of the moment equilibrium:
```math
- \\boldsymbol{\\sigma} \\cdot \\boldsymbol{\\nabla} = \\boldsymbol{0}
```

with boundary conditions on the boundary ``\\Gamma = \\Gamma_\\text{D} \\cup \\Gamma_\\text{N}``

```math
\\begin{aligned}
\\boldsymbol{u} &= \\boldsymbol{u}_\\text{p} \\quad \\text{on} \\quad \\Gamma_\\text{D} \\\\
\\boldsymbol{t} &= \\boldsymbol{t}_\\text{p} \\quad \\text{on} \\quad \\Gamma_\\text{N} \\\\
\\end{aligned}
```

where \$\\boldsymbol{u}\$ is the primary unknown field.
"""
struct Primal end

"""
    element_routine!(::Primal, ke, fe, cv, xe, material, thickness[, fv, grid, cellid, tₚ, faceset_name])

Compute the element stiffness matrix `ke` and the element load vector `fe`:

```math
\\begin{aligned}

\\left( ke \\right)_{ij} &= a\\left( \\boldsymbol{N}^{(\\rm u)}_i, \\, \\boldsymbol{N}^{(\\rm u)}_j \\right)
= \\int_\\Omega \\boldsymbol{\\nabla}^\\text{sym} \\boldsymbol{N}^{(\\rm u)}_i : \\boldsymbol{\\mathsf{E}} :
\\boldsymbol{\\nabla}^\\text{sym} \\boldsymbol{N}^{(\\rm u)}_j \\rm d \\Omega \\\\

\\left( fe \\right)_{i} &= l\\left( \\boldsymbol{N}^{(\\rm u)}_i \\right)
= \\int_{\\Gamma^{(\\rm u)}_\\text{N} } \\boldsymbol{N}^{(\\rm u)}_i \\cdot \\boldsymbol{t}_\\text{p} \\, \\rm d \\Gamma

\\end{aligned}
```

# Arguments:
- `ke`: element stiffness matrix
- `fe`: element force vector
- `cv`: CellVectorValues
- `xe`: element coordinate vector
- `material`: material parameters, for instant only `LinearElasticity` is allowed
- `thickness`: out-of-plane thickness

# Optional arguments for integrating the load vector:
- `fv`: FaceVectorValues
- `grid`: Grid
- `cellid`: global cell index of the current cell
- `tₚ`: traction vector
- `faceset_name`: name of the faceset over which the load vector should be integrated
"""
function element_routine!( 
    ::Primal,
    ke,
    fe,
    cv,
    xe,
    material,
    thickness,
    fv=nothing,
    grid = nothing,
    cellid = nothing,
    tₚ=nothing,
    faceset_name=nothing,
)
    fill!(ke, 0.0)
    fill!(fe, 0.0)

    # unpack variables
    E = material.Eᵉ
    
    reinit!(cv, xe)

    ndof = getnbasefunctions(cv)

    # assemble a(u, δu)
    for qp in 1:getnquadpoints(cv) # both cellvalues have the same number of qp
        detJ = getdetJdV(cv, qp)
        for j in 1:ndof
            ∇Nⱼᵘ_sym = shape_symmetric_gradient(cv, qp, j)
            for i in 1:ndof
                ∇Nᵢᵘ_sym = shape_symmetric_gradient(cv, qp, i)
                ke[i,j] += ∇Nᵢᵘ_sym ⊡ E ⊡ ∇Nⱼᵘ_sym * detJ * thickness
            end
        end
    end
    
    # assemble l(δu)
    if faceset_name !== nothing
        if cellid === nothing || grid === nothing || fv === nothing || tₚ === nothing
            error("You must give a cellid, grid, facevalues and a traction vector for integrating the load vector.")
        end
        for face in 1:nfaces(getcells(grid, cellid))
            if FaceIndex(cellid, face) ∈ getfaceset(grid, faceset_name)
                reinit!(fv, xe, face)
                for qp in 1:getnquadpoints(fv)
                    detJ = getdetJdV(fv, qp)
                    for i in 1:ndof
                        Nᵢᵘ = shape_value(fv, qp, i)
                        fe[i] += Nᵢᵘ ⋅ tₚ * detJ * thickness
                    end
                end
            end
        end
    end

    return ke, fe
end