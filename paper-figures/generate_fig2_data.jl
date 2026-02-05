using Pkg
Pkg.activate("..")
using Revise,CSV
using PythonPlot,Statistics,Distributions,DataFrames
PythonPlot.svg(true)
using SizeStructuredPopulations
using SpecialFunctions
include("./examples.jl")

models = [M1(), M2(),M0()]
init = vcat([0.0], [0.0, 0.0], zeros(3))
dt = 0.1
T_cell = 1000.5
T_lineage = 7000.0
T_population = 9.

function terminate(cell)
    return cell.time[end] > T_population
end




dfs = [] # each dataframe added to this should have columns, time,c,ensemble, model
for (i, m) in enumerate(models)
    print("Model $i: ")
    model = GrowthModel(m...)

    # lineage ----------------------------------------------------------------------
    cell = simulate_lineage(model, T_lineage, init, dt = 0.1);
    df = lineage_to_dataframe(cell);
    dfc_cell = DataFrame()
    dfc_cell.time = df.time
    lag = 10000
    dfc_cell.z1_x1 = [cor(df.z1[max(1,k-lag):k],df.x1[max(1,k-lag):k]) for k in 1:length(df.time)]
    dfc_cell.z1_x12 = [cor(df.z1[max(1,k-lag):k],df.x1[max(1,k-lag):k] .^2) for k in 1:length(df.time)]
    dfc_cell.z2_x12 = [cor(df.z2[max(1,k-lag):k],df.x1[max(1,k-lag):k] .^2) for k in 1:length(df.time)]
    dfc_cell.z0_x12 = [cor(df.z1[max(1,k-lag):k] .- df.z2[max(1,k-lag):k],df.x1[max(1,k-lag):k] .^2) for k in 1:length(df.time)]
    dfc_cell.model =  ones(length(dfc_cell.time))*i
    dfc_cell.ensemble .= "lineage"
    push!(dfs, dfc_cell)

    # population ----------------------------------------------------------------------
    n_samples = 40 # the number of initial cells 

    # we sample the initial cells from the lineage data whose distribution should be close to the population steady state
    init_data = hcat(zeros(n_samples),Matrix(df[sample(1:nrow(df), n_samples; replace=true), [:z1,:z2,:x1,:x2,:x3]]))
    init_cells = [generate_cell(model,init_data[i,:],T_cell) for i in 1:n_samples]
    [grow_tree!(r, model, terminate) for r in init_cells];

    # now make the population dataframe
    df_pops = []
    j = 1
    for r in init_cells
        df_pop = population_to_dataframe(r)
        df_pop.root_ind .= j
        j += 1
        push!(df_pops, df_pop)
    end
    df_pop = vcat(df_pops...);

    # now we get the correlations
    df_pop = df_pop[df_pop.time .< T_population,:]
    grouped = groupby(df_pop, :time)
    dfc_pop = combine(grouped, [:z1, :x1] => ((z, x) -> cor(z, x )) => :z1_x1,
        [:z1, :x1] => ((z, x) -> cor(z, x .^2 )) => :z1_x12,
        [:z2, :x1] => ((z, x) -> cor(z, x .^2 )) => :z2_x12,
        [:z2,:z1, :x1] => ((z2,z1, x) -> cor(z1 - z2, x .^2 )) => :z0_x12)
    dfc_pop.model = ones(length(dfc_pop.time))*i
    dfc_pop.ensemble .= "population"
    push!(dfs, dfc_pop)
  
end

df = vcat(dfs...);
CSV.write("fig2_correlation_data.csv", df)