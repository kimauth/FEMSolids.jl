# find the longest edge of a cell
function longest_edge(grid, cellid)

    max_length = 0.0
    faceid = 0

    cell = getcells(grid, cellid)
    for f in 1:nfaces(cell)
        # assumes linear elements
        node_idx1, node_idx2 = Ferrite.faces(cell)[f]
        l = norm(grid.nodes[node_idx1].x - grid.nodes[node_idx2].x)
        if l > max_length
            max_length = l
            faceid = f
        end
    end

    return FaceIndex(cellid, faceid)
end

# create a new node on the mid-point the the face given by faceindex
function new_node(grid, faceindex)
    cellid, faceid = faceindex
    cell = getcells(grid, cellid)
    facenodes = Ferrite.faces(cell)[faceid]
    X = mean([grid.nodes[i].x for i in facenodes]) # hard-coded linear elements
    return Node(X)
end

# add nodes to new_nodes Vector and to split_edges dict
function add_nodes!(new_nodes, split_edges, grid, face::FaceIndex, nodeid)
    faces = Ferrite.full_neighborhood(grid, face, true)
    if !haskey(split_edges, first(faces)) # this face hasn't been assigned a new node yet
        for f in faces[2:end]
            @assert !haskey(split_edges, f) # the neighboring face should not be stored either
        end
        push!(new_nodes, new_node(grid, first(faces)))
        nodeid += 1
        merge!(split_edges, Dict(faces .=> nodeid))
    end
    return nodeid
end

"""
    refine(grid::Grid, cells_to_split)

Refine the given `grid` by applying the [Rivara algorithm](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620200412) where `cells_to_split` is the list of cells to define (should contain cell numbers).
The `grid` must include a `Topology` for the refinement to work. Construct the first grid as follows, all grids that are returned by `refine` include the topology information:

```julia
grid = generate_grid(Triangle, (nel_x, nel_y), lower_left, upper_right; ; build_topology=true)
```

This implementation is restricted to linear 2D Triangles.

In order to use this function, you need to use the CA4 branch of Ferrite.jl. You can add it as follows:
```
pkg> add Ferrite#CA4
```

!!! warning
    `nodesets`, `facesets` and `cellsets` are currently lost after refinement and must be reconstructed. They can be added to a grid in the following manner:
    ```julia
    addfaceset!(mesh, "left", x -> x[1] ≈ 0.0)
    addfaceset!(mesh, "bottom", x -> x[2] ≈ 0.0)
    addfaceset!(mesh, "Γₗ", x -> (x[1] <= h && x[2] ≈ h))
    ```
"""
function refine(grid::Grid{2,Triangle}, cells_to_split)

    # compute new nodes
    nnodes = getnnodes(grid)
    nodeid = nnodes
    new_nodes = Node{2,Float64}[]
    split_edges = Dict{FaceIndex, Int}()

    old_cell_list = Int[]
    cell_list = cells_to_split
    while !(old_cell_list == cell_list)
        for cellid in cell_list
            face = longest_edge(grid, cellid)
            nodeid = add_nodes!(new_nodes, split_edges, grid, face, nodeid)
        end
        old_cell_list = copy(cell_list)
        cell_list = unique!([f[1] for f in keys(split_edges)])
    end


    new_cells = Triangle[]
    # construct new cells
    for cellid in cell_list
        longest_face = longest_edge(grid, cellid)
        # n_face = Ferrite.full_neighborhood(grid, longest_face)[1] # assumes only one neighboring face
        node_longest_face = split_edges[longest_face]

        cell = getcells(grid, cellid)
        for f in 1:nfaces(cell)
            if f == longest_face[2]
                continue
            end
            face = Ferrite.faces(cell)[f]
            node_f = get(split_edges, FaceIndex(cellid, f), nothing)
            if node_f === nothing # no node on the face
                c = Triangle((face[1], face[2], node_longest_face))
                push!(new_cells, c)
            else
                c1 = Triangle((face[1], node_f, node_longest_face))
                c2 = Triangle((node_f, face[2], node_longest_face))
                push!(new_cells, c1)
                push!(new_cells, c2)
            end
        end
    end

    # construct new grid
    cells = vcat(deleteat!(grid.cells, sort!(cell_list)), new_cells)
    nodes = vcat(grid.nodes, new_nodes)

    new_grid = Grid(cells, nodes; topology=Ferrite.ExclusiveTopology(cells)) # currently looses facesets etc.
    return new_grid
end

"""
    linear_to_quadratic(dh_lin::DofHandler{2,T,G}, dh_quad::DofHandler{2,T,G}, a_lin) where {T, G<:Grid{2,Triangle}}

Interpolate the linear solution `a_lin` (associated with `dh_lin`) to the degrees of freedom of the quadratic `DofHandler` `dh_quad`.
Return the interpolated values as a vector ordered according to the dofs in `dh_quad`. Restricted to vector valued fields.
"""
function linear_to_quadratic(dh_lin::DofHandler{2,T,G}, dh_quad::DofHandler{2,T,G}, a_lin) where {T, G<:Grid{2,Triangle}}
    dim = 2

    qr = QuadratureRule{dim,RefTetrahedron,Float64}([NaN for i=1:6], Vec.([(1.0, 0.0), (0.0, 1.0), (0.0, 0.0), (0.5, 0.5), (0.0, 0.5), (0.5, 0.0)]))
    ip = Lagrange{dim, RefTetrahedron, 1}()
    cv_nodes = CellVectorValues(qr, ip)


    a_quad = Vector{Float64}(undef, ndofs(dh_quad))

    for cell in CellIterator(dh_quad)
        dofs = celldofs(cell)
        dofs_lin = celldofs(dh_lin, cellid(cell))
        for (node, dof) in enumerate(1:dim:length(dofs))
            dof1 = dofs[dof]
            a_quad[dof1:(dof1+dim-1)] = function_value(cv_nodes, node, a_lin[dofs_lin])
        end
    end

    return a_quad
end
