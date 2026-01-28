# AI Agent Instructions for SizeStructuredPopulations

This Julia package implements stochastic simulations of size-structured cell populations with growth fluctuations, supporting the research paper "Size-structured populations with growth fluctuations: Feynman–Kac formula and decoupling".

## Architecture Overview

The codebase follows a biological simulation pattern with 3 core abstractions:

- **`Cell`**: Mutable struct storing time-series data (`time`, `z`, `x`) and binary tree structure (`daughterL`, `daughterR`) 
- **`GrowthModel`**: Immutable struct containing 4 functions defining population dynamics: `lambda` (growth), `beta` (division rate), `h` (division kernel), `L` (phenotype evolution)
- **Simulation modes**: Lineage tracking (single branch) vs population trees (full binary trees)

Key architectural insight: Cells maintain their complete temporal trajectory while forming tree structures through division events.

## Core Patterns

### Model Definition Pattern
Models are defined as tuples of functions returned by constructors (see `paper-figures/examples.jl`):
```julia
function M1()
    λ(x) = 1.0 + x[1]  # Growth rate depends on phenotype
    β(x, z, t) = λ(x) * A / B  # Division rate (complex hazard function)
    h(z, x) = ([z[1] - log(2), 0], x)  # Division kernel (size reset + inheritance)
    L(x, t, dt) = -γ*x*dt .+ sqrt(dt*D) .*rand(Normal(0,1), 3)  # OU phenotype dynamics
    return λ, β, h, L  # Return tuple, not GrowthModel struct
end
```

### State Vector Convention
All simulations use consistent state encoding: `[t, z1, z2, x1, x2, x3, ...]`
- `z1`: log-size, `z2`: log-size at division (reset to 0 after division)
- `x`: phenotypic state vector (typically 3D for examples)

### Data Export Pattern  
Two DataFrame export functions with different tree traversal:
- `lineage_to_dataframe(root)`: Single lineage (daughterL chain only)
- `population_to_dataframe(root)`: Full population tree (both daughters)

Both automatically expand `x` vectors into separate columns (`x1`, `x2`, `x3`).

## Development Workflow

### Environment Setup
```julia
using Pkg; Pkg.activate(".")  # Always activate project environment
using Revise  # Essential for interactive development
using SizeStructuredPopulations
```

### Typical Analysis Flow
1. Activate project environment in notebooks with `Pkg.activate("..")`  
2. Define model functions (see `examples.jl` patterns)
3. Set initial conditions: `init = vcat([0.0], [0.0, 0.0], zeros(3))`
4. Run simulation: `simulate_lineage(model, Tmax, init)` 
5. Convert to DataFrame for analysis: `lineage_to_dataframe(root)`
6. Use PythonPlot with `PythonPlot.svg(true)` for publication figures

### File Organization
- `src/`: Core library (`model.jl`, `lineage.jl`, `population.jl`, utilities)
- `example-notebooks/`: Learning examples with simple models
- `paper-figures/`: Specific analyses for publication (complex models in `examples.jl`)

## Critical Dependencies
- **PythonPlot**: Primary plotting backend (not Plots.jl)
- **Distributions**: Essential for stochastic processes in model functions  
- **SpecialFunctions**: Required for `erf()` in division rate calculations
- **Revise**: Critical for interactive development workflow

## Common Gotchas
- Tree traversal inconsistency: `population.jl` uses `daughter_L`/`daughter_R` but should be `daughterL`/`daughterR`
- Division kernel `h()` must return tuple of vectors: `([new_z], [new_x])`
- Time integration uses simple Euler scheme with fixed `dt=0.01` 
- Random seeds are set per-notebook for reproducible figure generation
- Model constructors return function tuples, not GrowthModel structs (inconsistent pattern)

Focus on maintaining the biological simulation semantics and consistent state vector conventions when extending functionality.