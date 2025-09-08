module SizeStructuredPopulations
using CSV,DataFrames
include("model.jl")
include("lineage.jl")
include("population.jl")
include("utils.jl")
include("utils_plotting.jl")

export GrowthModel, Cell, generate_cell, simulate, simulate_lineage, lineage_to_dataframe, grow_tree!
export make_branches,population_to_dataframe

end