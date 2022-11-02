
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
    element_routine!(::Primal, ke, fe, fe_external, cv, xe, material, thickness[, fv, neumann_bc::Neumann])

Compute the element stiffness matrix `ke`, the internal force vector `fe` and the external 
force vector `fe_external`:

```math
\\begin{aligned}

\\left( \\texttt{ke} \\right)_{ij} &= a\\left( \\boldsymbol{N}^{(\\rm u)}_i, \\, \\boldsymbol{N}^{(\\rm u)}_j \\right)
= \\int_\\Omega \\boldsymbol{\\nabla}^\\text{sym} \\boldsymbol{N}^{(\\rm u)}_i : \\boldsymbol{\\mathsf{E}} :
\\boldsymbol{\\nabla}^\\text{sym} \\boldsymbol{N}^{(\\rm u)}_j \\rm d \\Omega \\\\

\\left( \\texttt{fe} \\right)_{i} &= a\\left( \\boldsymbol{u}, \\, \\boldsymbol{N}^{(\\rm u)}_j \\right)
= \\int_\\Omega \\boldsymbol{\\nabla}^\\text{sym} \\boldsymbol{u} : \\boldsymbol{\\mathsf{E}} :
\\boldsymbol{\\nabla}^\\text{sym} \\boldsymbol{N}^{(\\rm u)}_j \\rm d \\Omega \\\\

\\left( \\texttt{fe\\_external} \\right)_{i} &= l\\left( \\boldsymbol{N}^{(\\rm u)}_i \\right)
= \\int_{\\Gamma^{(\\rm u)}_\\text{N} } \\boldsymbol{N}^{(\\rm u)}_i \\cdot \\boldsymbol{t}_\\text{p} \\, \\rm d \\Gamma,

\\end{aligned}
```
where \$\\boldsymbol{N}^{(\\rm u)}_i\$ is the shape function for the i-th degree of freedom,
\$ \\mathsf{E} \$ is the 4-th order elastic stiffness tensor, \$\\boldsymbol{u}\$ is the 
displacement and \$\\boldsymbol{t}_\\mathrm{p}\$ is the prescribed traction.

# Arguments:
- `ke`: element stiffness matrix
- `fe`: element force vector
- `cv`: CellVectorValues
- `xe`: element coordinate vector
- `material`: material parameters, for instant only `LinearElasticity` is allowed
- `thickness`: out-of-plane thickness

# Optional arguments for integrating the load vector:
- `fv`: FaceVectorValues
- `neumann_bc`: `Neumann` boundary condition
"""
function element_routine!( 
    ::Primal,
    ke,
    fe,
    fe_external,
    cv,
    xe,
    material,
    thickness,
    fv = nothing,
    neumann_bc = nothing,
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
    
    # assemble l(δu) (no body loads)
    if neumann_bc !== nothing
        # destructure named fields from struct into variables
        (; f, faceset, nfaces, time, current_cellid) = neumann_bc
        for face in 1:nfaces
            if FaceIndex(current_cellid[], face) ∈ faceset # check face in on Neumann boundary
                reinit!(fv, xe, face)
                for qp in 1:getnquadpoints(fv)
                    x = spatial_coordinate(fv, qp, xe)
                    tₚ = f(x, time[])
                    detJ = getdetJdV(fv, qp)
                    for i in 1:ndof
                        Nᵢᵘ = shape_value(fv, qp, i)
                        fe_external[i] += Nᵢᵘ ⋅ tₚ * detJ * thickness
                    end
                end
            end
        end
    end

    return ke, fe, fe_external
end

