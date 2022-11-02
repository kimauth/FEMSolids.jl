@testset "CST element" begin
    # test stiffness matrix
    dim = 2
    material = LinearElasticity{2}(K=1.25, G=0.4)

    ip = Lagrange{dim, RefTetrahedron, 1}()
    qr = QuadratureRule{dim, RefTetrahedron}(1)
    cv = CellVectorValues(qr, ip)

    xe = [Vec(6., -4.), Vec(5., 3.), Vec(-4., -1.)]

    thickness = 2.0

    ke = Matrix{Float64}(undef, 6, 6)
    fe = Vector{Float64}(undef, 6)

    qr_face = QuadratureRule{dim-1, RefTetrahedron}(1)
    fv = FaceVectorValues(qr_face, ip)

    nodes = Node.(xe)
    cells = [Triangle((1,2,3))]
    grid = Grid(cells, nodes)
    addfaceset!(grid, "Γ", Set((FaceIndex(1,1),)))

    weak_form = Primal(material, thickness)
    tₚ = Vec(1.0, 1.0)

    neumann_bc = Neumann((x,t)->tₚ, getfaceset(grid, "Γ"), nfaces(getcells(grid, 1)))
    update!(neumann_bc)
    update_cell!(neumann_bc, 1)
    
    ue = zeros(6)

    element_routine!(weak_form, ke, fe, fe, cv, xe, ue, fv, neumann_bc)

    ke_reference = [
        0.9095 -0.7433 -0.2179 0.4259 -0.6915 0.3174
        -0.7433 2.2515 -0.1575 -2.3239 0.9007 0.0724
        -0.2179 -0.1575 0.8366 0.6194 -0.6187 -0.4619
        0.4259 -2.3239 0.6194 2.7154 -1.0453 -0.3915
        -0.6915 0.9007 -0.6187 -1.0453 1.3102 0.1445
        0.3174 0.0724 -0.4619 -0.3915 0.1445 0.3192
        ]
    @test round.(ke, digits=4) ≈ ke_reference

    fe_reference = tₚ*thickness*norm(xe[2]-xe[1]) / 2
    @test fe[1:4] ≈ vcat(fe_reference, fe_reference)
end

@testset "Neumann boundary condition" begin

    grid = generate_grid(Quadrilateral, (1,1))

    xe = [getnodes(grid, n).x for n in getcells(grid, 1).nodes]

    ip = Lagrange{2,RefCube,1}()
    qr = QuadratureRule{2,RefCube}(2)
    cv = CellVectorValues(qr, ip)

    qr_face = QuadratureRule{1,RefCube}(1)
    fv = FaceVectorValues(qr_face, ip)

    ke = zeros(8,8)
    fe = zeros(8)

    material = LinearElasticity{2}(G=1.0, K=1.0)
    weak_form = Primal(material, 1.0)

    tₚ = Vec(1.0, 2.0)
    
    neumann_bc = Neumann((x,t)->tₚ, getfaceset(grid, "top"), nfaces(getcells(grid, 1)))
    update!(neumann_bc)
    update_cell!(neumann_bc, 1)
    
    ue = zeros(8)

    element_routine!(weak_form, ke, fe, fe, cv, xe, ue, fv, neumann_bc)

    @test fe == [0., 0., 0., 0., 1., 2., 1., 2.]
end

   
