# find the longest edge of a cell
function longest_edge(grid, cellid)

    max_length = 0.0
    edgeid = 0

    cell = getcells(grid, cellid)
    for f in 1:nfacets(cell)
        # assumes linear 2d elements (facet = edge)
        node_idx1, node_idx2 = Ferrite.facets(cell)[f]
        l = norm(grid.nodes[node_idx1].x - grid.nodes[node_idx2].x)
        if l > max_length
            max_length = l
            edgeid = f
        end
    end

    return EdgeIndex(cellid, edgeid)
end

# create a new node on the mid-point the the facet given by facetindex
function new_node(grid, edgeindex)
    cellid, edgeid = edgeindex
    cell = getcells(grid, cellid)
    edgenodes = Ferrite.facets(cell)[edgeid]
    x = mean(grid.nodes[i].x for i in edgenodes) # hard-coded linear elements
    return Node(x)
end

# add nodes to new_nodes Vector and to split_edges dict
function add_nodes!(new_nodes, split_edges, grid, topology, edge::EdgeIndex, nodeid)
    edges = getneighborhood(topology, grid, edge, true)
    if !haskey(split_edges, first(edges)) # this edge hasn't been assigned a new node yet
        for e in edges[2:end]
            @assert !haskey(split_edges, e) # the neighboring facet should not be stored either
        end
        push!(new_nodes, new_node(grid, first(edges)))
        nodeid += 1
        merge!(split_edges, Dict(edges .=> nodeid))
    end
    return nodeid
end

"""
    refine(grid::Grid, cells_to_split)

Refine the given `grid` by applying the [Rivara algorithm](https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.1620200412) where `cells_to_split` is the list of cells to define (should contain cell numbers).
This implementation is restricted to linear 2D Triangles.

!!! warning
    `nodesets`, `facesets` and `cellsets` are currently lost after refinement and must be reconstructed. They can be added to a grid in the following manner:
    ```julia
    addfaceset!(mesh, "left", x -> x[1] ≈ 0.0)
    addfaceset!(mesh, "bottom", x -> x[2] ≈ 0.0)
    addfaceset!(mesh, "Γₗ", x -> (x[1] <= h && x[2] ≈ h))
    ```
"""
function refine(grid::Grid{2,Triangle}, cells_to_split=1:getncells(grid))
    # compute topology
    topology = ExclusiveTopology(grid)

    # Should be possible to work with the original grid directly,
    # but to be on the safe side to not alter the original implementation
    tmp_grid = Grid(copy(grid.cells), copy(grid.nodes))

    # compute new nodes
    nnodes = getnnodes(tmp_grid)
    nodeid = nnodes
    new_nodes = Node{2,Float64}[]
    split_edges = Dict{EdgeIndex, Int}()

    old_cell_list = Int[]
    cell_list = cells_to_split
    while !(old_cell_list == cell_list)
        for cellid in cell_list
            edge = longest_edge(tmp_grid, cellid)
            nodeid = add_nodes!(new_nodes, split_edges, tmp_grid, topology, edge, nodeid)
        end
        old_cell_list = copy(cell_list)
        cell_list = unique!([f[1] for f in keys(split_edges)])
    end


    new_cells = Triangle[]
    # construct new cells
    for cellid in cell_list
        the_longest_edge = longest_edge(tmp_grid, cellid)
        node_longest_edge = split_edges[the_longest_edge]

        cell = getcells(tmp_grid, cellid)
        for f in 1:nfacets(cell)
            if f == the_longest_edge[2]
                continue
            end
            facet = Ferrite.facets(cell)[f]
            node_f = get(split_edges, EdgeIndex(cellid, f), nothing)
            if node_f === nothing # no node on the facet
                c = Triangle((facet[1], facet[2], node_longest_edge))
                push!(new_cells, c)
            else
                c1 = Triangle((facet[1], node_f, node_longest_edge))
                c2 = Triangle((node_f, facet[2], node_longest_edge))
                push!(new_cells, c1)
                push!(new_cells, c2)
            end
        end
    end

    # construct new grid
    cells = vcat(deleteat!(tmp_grid.cells, sort!(cell_list)), new_cells)
    nodes = vcat(tmp_grid.nodes, new_nodes)

    new_grid = Grid(cells, nodes)
    return new_grid
end

"""
    transfer_solution(dh_from::DofHandler, dh_to::DofHandler, a_from::AbstractVector)

Interpolate the solution `a_from` (associated with `dh_from`) to the degrees of freedom of the `DofHandler` `dh_to`.
Return the interpolated values as a vector ordered according to the dofs in `dh_to`.
"""
function transfer_solution(dh_from::DofHandler, dh_to::DofHandler, a_from::AbstractVector)
    # Checks
    if Ferrite.get_grid(dh_from) !== Ferrite.get_grid(dh_to)
        error("Both dofhandlers should have the same underlying grid")
    end
    if any(key ∉ Ferrite.getfieldnames(dh_from) for key in Ferrite.getfieldnames(dh_to))
        error("All fields in `dh_to` must be in `dh_from`")
    end
    @assert length(a_from) == ndofs(dh_from)
    
    # Would be possible to lift by generalizing the code
    if length(dh_from.subdofhandlers) != length(dh_to.subdofhandlers) != 1
        error("Both dh_from and dh_to can only have a single subdofhandler")
    end
    sdh_from = only(dh_from.subdofhandlers)
    sdh_to = only(dh_to.subdofhandlers)
    
    ip_geo = geometric_interpolation(getcelltype(sdh_to))
    a_to = similar(a_from, ndofs(dh_to))
    for fieldname in Ferrite.getfieldnames(sdh_to)
        ip_from = Ferrite.getfieldinterpolation(sdh_from, fieldname)
        ip_to = Ferrite.getfieldinterpolation(sdh_to, fieldname)
        ξ_to = Ferrite.reference_coordinates(Ferrite.get_base_interpolation(ip_to))
        qr = QuadratureRule{getrefshape(ip_to)}(fill(NaN, length(ξ_to)), ξ_to)
        cv_nodes = CellValues(qr, ip_from, ip_geo; update_gradients = false, update_detJdV = false)
        _transfer_solution!(a_to, sdh_to, a_from, sdh_from, cv_nodes, dof_range(sdh_to, fieldname), dof_range(sdh_from, fieldname))
    end
    return a_to
end

function _transfer_solution!(a_to::AbstractVector{T}, sdh_to, a_from, sdh_from, cv, dr_to, dr_from) where {T}
    ae_from = zeros(T, length(dr_from))
    cell_from = CellCache(sdh_from)
    for cell_to in CellIterator(sdh_to)
        reinit!(cell_from, cellid(cell_to))
        reinit!(cv, cell_to)
        dofs_to = celldofs(cell_to)
        copy!(ae_from, view(a_from, celldofs(cell_from)))
        i = 0
        for dof_location in 1:getnquadpoints(cv)
            u = function_value(cv, dof_location, ae_from, dr_from)
            r = dr_to[i .+ (1:length(u))]
            for (d, uval) in zip(r, u)
                a_to[dofs_to[d]] = uval
            end
            i += length(u)
        end
    end
    return a_to
end
