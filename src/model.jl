
mutable struct Cell
    time::Vector{Float64}
    z::Vector{Vector{Float64}}
    x::Vector{Vector{Float64}}
    daughterL::Union{Cell, Nothing}
    daughterR::Union{Cell, Nothing}
end


# Outer constructor: no daughters specified
Cell(time, z, x) = Cell(time, z, x, nothing, nothing)
Cell(time, z, x, daughterL, daughterR) = Cell(time, z, x, daughterL, daughterR)


struct GrowthModel
    lambda::Function           # Growth rate function: λ(x, t)
    beta::Function             # Division rate function: β(x,z, t)
    h::Function                # Division kernel: h(z′ | z)
    L::Function                # phenotype evolution function: L(x, t, dt)
end




function generate_cell(model::GrowthModel, init::Vector{Float64},Tmax::Float64; dt=0.01)
    t = init[1]
    z0 = init[2:3]
    x0 = init[4:end]

    times = [t]
    z = [copy(z0)]
    x = [copy(x0)]


    div = false
    while (!div) && (t < Tmax)
        t += dt
        λ = model.lambda(x[end])
        z_new = z[end] .+  λ .* dt
        x_new = x[end] .+  model.L(x[end], t, dt)
        #x_new = [max(0.01,xi) for xi in xnew]

        push!(times, t)
        push!(z, z_new)
        push!(x, x_new)

        div_prob = model.beta(x[end], z[end], t) * dt
        if rand() < div_prob
            div = true
        end
    end

    return Cell(times, z, x, nothing, nothing)  # Using outer constructor which handles daughter cells
end


function simulate(model::GrowthModel, Tmax::Float64, ensemble::Symbol, init::Vector{Float64})
    if ensemble == :lineage
        return simulate_lineage(model, Tmax, init)
    else
        return simulate_lineage(model, Tmax, init)
    end
end