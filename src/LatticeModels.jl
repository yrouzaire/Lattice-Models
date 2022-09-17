include("lattices.jl")
include("models.jl")
include("init_visu.jl")

include("core_methods.jl")
include("defects_methods.jl")

include("misc.jl")
include("measurements.jl")

using BenchmarkTools

# using Flux:onecold, Chain, Dense, softmax
# global const NN = load("NN_all_12_defects.jld2","NN")
# global const possible_defects = load("NN_all_12_defects.jld2","possible_defects")
