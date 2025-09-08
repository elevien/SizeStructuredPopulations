function make_branches(root::Cell,dt,max_depth)
    branches = Vector{Tuple{Float64, Float64, Float64, Float64}}()

    function traverse(cell::Cell, x0::Float64, y0::Float64, depth::Int)
        # Compute vertical branch (lifespan)
        dt_life = (cell.time[end] - cell.time[1])
        y1 = y0 - dt_life
        push!(branches, (x0, y0, x0, y1))
        
        dx = 2.0^(-depth)

        # Horizontal links to daughters
        y_daughter = y1  # All daughters start at y1
        if depth <= max_depth
            # Check for daughters
            if cell.daughterL !== nothing 
                # DaughterL
                xL = x0 - dx
                push!(branches, (x0, y1, xL, y_daughter))
                traverse(cell.daughterL, xL, y_daughter, depth + 1)
            end
            if cell.daughterR !== nothing
                # DaughterR
                xR = x0 + dx
                push!(branches, (x0, y1, xR, y_daughter))
                traverse(cell.daughterR, xR, y_daughter, depth + 1)
            end
        end
    end

    traverse(root, 0.0, 0.0, 2)
    return branches
end