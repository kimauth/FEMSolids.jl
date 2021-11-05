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

    tₚ = Vec(1.0, 1.0)

    element_routine!(Primal(), ke, fe, cv, xe, material, thickness, fv, grid, 1, tₚ, "Γ")

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