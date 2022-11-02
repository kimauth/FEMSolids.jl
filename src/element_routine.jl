
"""
    Primal(material, thickness::Float64)

Represent the weak form of the moment equilibrium in primal format:
```math
\\int_\\Omega \\boldsymbol{\\sigma} : \\boldsymbol{\\varepsilon}(\\delta \\boldsymbol{u}) \\, \\mathrm{d}\\Omega =
\\int_\\Omega \\boldsymbol{b} \\cdot \\delta\\boldsymbol{u} \\, \\mathrm{d}\\Omega + \\int_{\\partial\\Omega_\\mathrm{N}} \\boldsymbol{t}_\\mathrm{p} \\cdot \\delta\\boldsymbol{u} \\, \\mathrm{d}\\Gamma,
```
where \$\\boldsymbol{\\sigma}\$ is the stress tensor, \$ \\boldsymbol{\\varepsilon} \$ is 
the strain tensor, \$ \\boldsymbol{b} \$ is the vector of body forces, \$ \\boldsymbol{t}_\\mathrm{p} \$
is the prescribed boundary traction and \$ \\delta \\boldsymbol{u} \$ is the test function.

Boundary conditions on the Dirichlet boundary ``\\Gamma_\\text{D}`` and the Neumann
boundary ``\\Gamma_\\text{N}`` are respectively given as:

```math
\\begin{aligned}
\\boldsymbol{u} &= \\boldsymbol{u}_\\text{p} \\quad \\text{on} \\quad \\Gamma_\\text{D} \\\\
\\boldsymbol{\\sigma} \\cdot \\boldsymbol{n} &= \\boldsymbol{t}_\\text{p} \\quad \\text{on} \\quad \\Gamma_\\text{N} \\\\
\\end{aligned}
```
where \$\\boldsymbol{u}\$ is the displacement field.

# Arguments:
- `materal`: Material model for the stress-strain relation
- `thickness`: Cross-sectional area / out-of-plane thickness for 1D / 2D problems 
"""
struct Primal{M}
    material::M
    thickness::Float64
end


"""
    element_routine!(weak_form::Primal, ke, fe, fe_external, cv, xe, ue [, fv, neumann_bc::Neumann])

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
- `weak_form`: weak form representation
- `ke`: element stiffness matrix
- `fe`: element force vector
- `cv`: CellVectorValues
- `xe`: element coordinate vector
- `ue`: element displacement vector

# Optional arguments for integrating the load vector:
- `fv`: FaceVectorValues
- `neumann_bc`: `Neumann` boundary condition
"""
function element_routine!( 
    weak_form::Primal,
    ke,
    fe,
    fe_external,
    cv,
    xe,
    ue,
    fv = nothing,
    neumann_bc = nothing,
)
    fill!(ke, 0.0)
    fill!(fe, 0.0)
    fill!(fe_external, 0.0)

    # unpack variables
    (; material, thickness) = weak_form
    E = material.Eᵉ
    
    reinit!(cv, xe)

    ndof = getnbasefunctions(cv)

    # assemble a(u, δu)
    for qp in 1:getnquadpoints(cv) # both cellvalues have the same number of qp
        detJ = getdetJdV(cv, qp)
        ε = function_symmetric_gradient(cv, qp, ue)
        for j in 1:ndof
            ∇Nⱼᵘ_sym = shape_symmetric_gradient(cv, qp, j)
            fe[j] += ε ⊡ E ⊡ ∇Nⱼᵘ_sym
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

