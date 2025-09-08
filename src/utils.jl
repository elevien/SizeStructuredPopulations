function lineage_to_dataframe(root::Cell)
    rows = []

    cell = root
    ind = 1.
    while cell !== nothing
        n = length(cell.time)
        d = length(cell.x[1])  # dimension of phenotype vector
        for i in 1:n
            row = (
                time = cell.time[i],
                z1 = cell.z[i][1],
                z2 = cell.z[i][2],
                x = cell.x[i],
                position = ind
            )
            push!(rows, row)
        end
        ind = ind + 1 
        cell = cell.daughterL
    end

    # Convert to DataFrame and expand `x` vector into separate columns
    df = DataFrame(rows)
    d = length(df.x[1])
    for j in 1:d
        df[!, "x$j"] = getindex.(df.x, j)
    end
    select!(df, Not(:x))  # remove the original `x` column
    return df
end

function lineage_to_matrix(root::Cell)
    df = lineage_to_dataframe(root)
    return Matrix(df)
end



function population_to_dataframe(root::Cell)
    rows = []

    function traverse(cell::Cell, branch::String)
        n = length(cell.time)
        d = length(cell.x[1])
        for i in 1:n
            row = (
                time = cell.time[i],
                z1 = cell.z[i][1],
                z2 = cell.z[i][2],
                x = cell.x[i],
                branch = branch,
            )
            push!(rows, row)
        end
        if cell.daughterL !== nothing
            traverse(cell.daughterL, branch * "0")
        end
        if cell.daughterR !== nothing
            traverse(cell.daughterR, branch * "1")
        end
    end

    traverse(root, "")
    df = DataFrame(rows)

    # Expand `x` vector into individual columns
    d = length(df.x[1])
    for j in 1:d
        df[!, "x$j"] = getindex.(df.x, j)
    end
    select!(df, Not(:x))  # remove the original `x` column

    return df
end