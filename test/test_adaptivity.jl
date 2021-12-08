@testset "grid refinement" begin
    #### test grid refinement against Figure 6 from https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620200412

    # pvd = paraview_collection("test")
    # test with grid from paper
    nodes = Node.([(0.0, 0.0), (2., 1.), (-1.5, 2.2), (2.4, 2.1), (2.3, 4.0), (0.3, 4.2), (-1.1, 4.1)])
    cells = Triangle.([(1,2,3), (2,4,3), (3,4,6), (4,5,6), (3,6,7)])
    grid = Grid(cells, nodes; topology=Ferrite.ExclusiveTopology(cells))
    # pvd[0] = vtk_grid("initial_grid", grid)

    cells_to_split = [5]
    grid_1 = refine(grid, cells_to_split)
    @test grid_1.cells == [  Triangle((1, 2, 3))
                            Triangle((4, 5, 6))
                            Triangle((4, 6, 9))
                            Triangle((6, 8, 9))
                            Triangle((8, 3, 9))
                            Triangle((2, 4, 9))
                            Triangle((3, 2, 9))
                            Triangle((6, 7, 8))
                            Triangle((7, 3, 8))]
    # pvd[1] = vtk_grid("1st_grid", grid_1)

    cells_to_split = [8,9]
    grid_2 = refine(grid_1, cells_to_split)
    @test grid_2.cells == [ Triangle((1, 2, 3))
                            Triangle((4, 5, 6))
                            Triangle((4, 6, 9))
                            Triangle((6, 8, 9))
                            Triangle((8, 3, 9))
                            Triangle((2, 4, 9))
                            Triangle((3, 2, 9))
                            Triangle((7, 8, 10))
                            Triangle((8, 6, 10))
                            Triangle((3, 8, 11))
                            Triangle((8, 7, 11))]
    # pvd[2] = vtk_grid("2nd_grid", grid_2)

    cells_to_split = [8,9,10,11]
    grid_3 = refine(grid_2, cells_to_split)
    @test grid_3.cells == [ Triangle((2, 4, 9))
                            Triangle((2, 9, 17))
                            Triangle((9, 16, 17))
                            Triangle((16, 3, 17))
                            Triangle((6, 10, 13))
                            Triangle((10, 8, 13))
                            Triangle((4, 5, 18))
                            Triangle((5, 6, 18)) 
                            Triangle((8, 10, 12))
                            Triangle((10, 7, 12))
                            Triangle((1, 2, 17)) 
                            Triangle((3, 1, 17)) 
                            Triangle((6, 13, 15))
                            Triangle((13, 8, 15))
                            Triangle((8, 9, 15)) 
                            Triangle((7, 11, 12))
                            Triangle((11, 8, 12))
                            Triangle((8, 11, 14))
                            Triangle((11, 3, 14))
                            Triangle((8, 14, 16))
                            Triangle((14, 3, 16))
                            Triangle((9, 8, 16)) 
                            Triangle((6, 15, 18))
                            Triangle((15, 9, 18))
                            Triangle((9, 4, 18))]
    # pvd[3] = vtk_grid("3rd_grid", grid_3)

    # vtk_save(pvd)
end

@testset "linear to quadratic interpolation" begin

    grid = generate_grid(Triangle, (1,1))

    dh_lin = DofHandler(grid)
    ip_u = Lagrange{2,RefTetrahedron,1}()
    push!(dh_lin, :u, 2, ip_u)
    close!(dh_lin)

    dh_quad = DofHandler(grid)
    ip_u = Lagrange{2,RefTetrahedron,2}()
    push!(dh_quad, :u, 2, ip_u)
    close!(dh_quad)

    a_lin = [1., 1., 2., 2., 2., 2., 4., 4.]
    a_quad = linear_to_quadratic(dh_lin, dh_quad, a_lin)

    @test a_quad == [1., 1., 2., 2., 2., 2., 1.5, 1.5, 2., 2., 1.5, 1.5, 4., 4., 3., 3., 3., 3.]
end