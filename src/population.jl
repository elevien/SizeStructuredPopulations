
"""
    grow_tree!(node, model, terminate; dt=0.01, parallel_depth=0, depth=0)

Grow a population tree by recursively generating daughter cells until `terminate` is true.
Set `parallel_depth` to the number of initial tree levels to spawn in parallel (0 keeps serial
execution). The `depth` keyword tracks recursion depth internally and should be left at its
default when calling.
"""
function grow_tree!(node,model::GrowthModel,terminate;dt=0.01,parallel_depth=0,depth=0)
    parallel_depth < 0 && throw(ArgumentError("parallel_depth must be >= 0"))
    if !terminate(node)
        z_parent = node.z[end]
        x_parent = node.x[end]
        t = node.time[end] + dt
    
        # Get daughter state from division kernel
        z_new, x_new = model.h(z_parent, x_parent)
        init = vcat(t, z_new, x_new)
        node.daughterL = generate_cell(model, init, Inf,dt=dt)
        node.daughterR = generate_cell(model, init, Inf,dt=dt)
        if depth < parallel_depth
            Threads.@sync begin
                Threads.@spawn grow_tree!(node.daughterL,model,terminate,dt=dt,parallel_depth=parallel_depth,depth=depth + 1)
                Threads.@spawn grow_tree!(node.daughterR,model,terminate,dt=dt,parallel_depth=parallel_depth,depth=depth + 1)
            end
        else
            grow_tree!(node.daughterL,model,terminate,dt=dt,parallel_depth=parallel_depth,depth=depth + 1)
            grow_tree!(node.daughterR,model,terminate,dt=dt,parallel_depth=parallel_depth,depth=depth + 1)
        end
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
