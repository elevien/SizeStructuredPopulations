function simulate_population(model::GrowthModel, Td::Float64, Tmax::Float64, Nmax::Int, init::Vector{Float64}; dt=0.01)
    # Parameters to pass down
    params = (model=model, dt=dt)

    # Termination condition: stop each tree at t >= Td
    terminate = cell -> cell.time[end] >= Td

    # Initial population of Nmax root cells
    population = [generate_cell(model, init, Tmax, dt=dt) for _ in 1:Nmax]

    # Grow each tree to time Td
    grow_forest!(population, terminate, params, (model, init, Tmax, dt) -> generate_cell(model, init, Tmax, dt=dt))

    # Sample new roots from leaf pool
    leaves = []
    for cell in population
        collect_leaves!(leaves, cell)
    end
    new_roots = rand(leaves, min(Nmax, length(leaves)))

    # Loop through dilution cycles until Tmax
    current_time = Td
    while current_time < Tmax
        terminate = cell -> cell.time[end] >= current_time + Td
        grow_forest!(new_roots, terminate, params, (model, init, Tmax, dt) -> generate_cell(model, init, Tmax, dt=dt))

        leaves = []
        for cell in new_roots
            collect_leaves!(leaves, cell)
        end
        new_roots = rand(leaves, min(Nmax, length(leaves)))

        current_time += Td
    end

    return new_roots  # forest of root cells
end

function grow_forest!(nodes::Vector{Cell}, terminate, params, generator)
    for node in nodes
        grow_tree!(node, terminate, params, generator)
    end
end

function grow_tree!(node::Cell, terminate, params, generator)
    if terminate(node)
        return
    end

    z_parent = node.z[end]
    x_parent = node.x[end]
    t = node.time[end]

    z1, x1 = params.model.h(z_parent, x_parent)
    init1 = vcat(t, z1, x1)
    node.daughterL = generator(params.model, init1, Tmax, dt=params.dt)
    grow_tree!(node.daughterL, terminate, params, generator)

    z2, x2 = params.model.h(z_parent, x_parent)
    init2 = vcat(t, z2, x2)
    node.daughterR = generator(params.model, init2, Tmax, dt=params.dt)
    grow_tree!(node.daughterR, terminate, params, generator)
end

function collect_leaves!(leaves::Vector{Cell}, node::Cell)
    if node.daughter_L === nothing && node.daughter_R === nothing
        push!(leaves, node)
    else
        if node.daughter_L !== nothing
            collect_leaves!(leaves, node.daughter_L)
        end
        if node.daughter_R !== nothing
            collect_leaves!(leaves, node.daughter_R)
        end
    end
end