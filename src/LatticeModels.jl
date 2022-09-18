include("lattices.jl")
include("models.jl")
include("init_visu.jl")

include("core_methods.jl")
include("defects_methods.jl")

include("misc.jl")
include("measurements.jl")

using BenchmarkTools

global const WINDOW = 7
# using Flux:onecold, Chain, Dense, softmax
# global const NN_positive = load("NN_positive_12_defects_N1000_W7.jld2","NN")
# global const possible_positive_defects = load("NN_positive_12_defects_N1000_W7.jld2","possible_defects")
# global const NN_negative = load("NN_negative_12_defects_N1000_W7.jld2","NN")
# global const possible_negative_defects = load("NN_negative_12_defects_N1000_W7.jld2","possible_defects")
