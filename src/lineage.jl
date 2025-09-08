

function simulate_lineage(model::GrowthModel, Tmax::Float64, init::Vector{Float64}; dt=0.01)
    # Generate the root cell
    root = generate_cell(model, init, Tmax, dt=dt)
    current = root

    while current.time[end] < Tmax
        # Get final state from current cell
        z_parent = current.z[end]
        x_parent = current.x[end]
        t = current.time[end] + dt

        # Get daughter state from division kernel
        z_new, x_new = model.h(z_parent, x_parent)
        init = vcat(t, z_new, x_new)

        # Generate next cell and attach as daughter_L
        daughter = generate_cell(model, init, Tmax,dt=dt)
        current.daughterL = daughter
        current = daughter
    end

    return root
end


