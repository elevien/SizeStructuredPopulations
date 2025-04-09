module SizeStructuredPopulations
using CSV,DataFrames
include("model.jl")
include("lineage.jl")
include("population.jl")
include("utils.jl")

export GrowthModel, Cell, generate_cell, simulate, simulate_lineage, lineage_to_dataframe, simulate_population

end