var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = FEMSolids","category":"page"},{"location":"#FEMSolids","page":"Home","title":"FEMSolids","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for FEMSolids.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package provides routines that are needed for solving the assignments in FEM Solids.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [FEMSolids]","category":"page"},{"location":"#FEMSolids.LinearElasticity-Union{Tuple{}, Tuple{dim}} where dim","page":"Home","title":"FEMSolids.LinearElasticity","text":"LinearElasticity{dim}(;G, K) where dim\n\nConstruct a LinearElasticity material based on the shear modulus G and the bulk modulus K.\n\n\n\n\n\n","category":"method"},{"location":"#FEMSolids.Primal","page":"Home","title":"FEMSolids.Primal","text":"Primal\n\nThis struct is used for dispatching which element_routine! is used when you add your own element routines. Primal refers to the primal format of the moment equilibrium:\n\n- boldsymbolsigma cdot boldsymbolnabla = boldsymbol0\n\nwith boundary conditions on the boundary Gamma = Gamma_textD cup Gamma_textN\n\nbeginaligned\nboldsymbolu = boldsymbolu_textp quad texton quad Gamma_textD \nboldsymbolt = boldsymbolt_textp quad texton quad Gamma_textN \nendaligned\n\nwhere boldsymbolu is the primary unknown field.\n\n\n\n\n\n","category":"type"},{"location":"#FEMSolids.element_routine!","page":"Home","title":"FEMSolids.element_routine!","text":"element_routine!(::Primal, ke, fe, cv, xe, material, thickness[, fv, grid, cellid, tₚ, faceset_name])\n\nCompute the element stiffness matrix ke and the element external load vector fe:\n\nbeginaligned\n\nleft( ke right)_ij = aleft( boldsymbolN^(rm u)_i  boldsymbolN^(rm u)_j right)\n= int_Omega boldsymbolnabla^textsym boldsymbolN^(rm u)_i  boldsymbolmathsfE \nboldsymbolnabla^textsym boldsymbolN^(rm u)_j rm d Omega \n\nleft( fe right)_i = lleft( boldsymbolN^(rm u)_i right)\n= int_Gamma^(rm u)_textN  boldsymbolN^(rm u)_i cdot boldsymbolt_textp  rm d Gamma\n\nendaligned\n\nArguments:\n\nke: element stiffness matrix\nfe: element force vector\ncv: CellVectorValues\nxe: element coordinate vector\nmaterial: material parameters, for instant only LinearElasticity is allowed\nthickness: out-of-plane thickness\n\nOptional arguments for integrating the load vector:\n\nfv: FaceVectorValues\ngrid: Grid\ncellid: global cell index of the current cell\ntₚ: traction vector\nfaceset_name: name of the faceset over which the load vector should be integrated\n\n\n\n\n\n","category":"function"},{"location":"#FEMSolids.linear_to_quadratic-Tuple{DofHandler{2, Triangle}, DofHandler{2, Triangle}, Any}","page":"Home","title":"FEMSolids.linear_to_quadratic","text":"linear_to_quadratic(dh_lin::DofHandler{2,Triangle}, dh_quad::DofHandler{2,Triangle}, a_lin)\n\nInterpolate the linear solution a_lin (associated with dh_lin) to the degrees of freedom of the quadratic DofHandler dh_quad. Return the interpolated values as a vector ordered according to the dofs in dh_quad. Restricted to vector valued fields.\n\n\n\n\n\n","category":"method"},{"location":"#FEMSolids.refine-Tuple{Grid{2, Triangle}, Any}","page":"Home","title":"FEMSolids.refine","text":"refine(grid::Grid, cells_to_split)\n\nRefine the given grid by applying the Rivara algorithm where cells_to_split is the list of cells to define (should contain cell numbers). The grid must include a Topology for the refinement to work. Construct the first grid as follows, all grids that are returned by refine include the topology information:\n\ngrid = generate_grid(Triangle, (nel_x, nel_y), lower_left, upper_right; ; build_topology=true)\n\nThis implementation is restricted to linear 2D Triangles.\n\nIn order to use this function, you need to use the CA4 branch of Ferrite.jl. You can add it as follows:\n\npkg> add Ferrite#CA4\n\nwarning: Warning\nnodesets, facesets and cellsets are currently lost after refinement and must be reconstructed. They can be added to a grid in the following manner:addfaceset!(mesh, \"left\", x -> x[1] ≈ 0.0)\naddfaceset!(mesh, \"bottom\", x -> x[2] ≈ 0.0)\naddfaceset!(mesh, \"Γₗ\", x -> (x[1] <= h && x[2] ≈ h))\n\n\n\n\n\n","category":"method"}]
}
