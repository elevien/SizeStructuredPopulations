
function grow_tree!(node,model::GrowthModel,terminate;dt=0.01)
    if !terminate(node)
        z_parent = node.z[end]
        x_parent = node.x[end]
        t = node.time[end] + dt
    
        # Get daughter state from division kernel
        z_new, x_new = model.h(z_parent, x_parent)
        init = vcat(t, z_new, x_new)
        node.daughterL = generate_cell(model, init, Inf,dt=dt)
        node.daughterR = generate_cell(model, init, Inf,dt=dt)
        grow_tree!(node.daughterL,model,terminate,dt=dt)
        grow_tree!(node.daughterR,model,terminate,dt=dt)
    end
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