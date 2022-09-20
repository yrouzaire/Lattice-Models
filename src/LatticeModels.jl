include("lattices.jl")
include("models.jl")
include("init_visu.jl")

include("core_methods.jl")
include("defects_methods.jl")

include("misc.jl")
include("measurements.jl")

using BenchmarkTools

# global const WINDOW = 7
# using Flux:onecold, Chain, Dense, softmax
# global const NN_positive12 = load("NN_positive_12_defects_N1000_W7.jld2","NN")
# global const possible_positive12_defects = load("NN_positive_12_defects_N1000_W7.jld2","possible_defects")
# global const NN_negative12 = load("NN_negative_12_defects_N1000_W7.jld2","NN")
# global const possible_negative12_defects = load("NN_negative_12_defects_N1000_W7.jld2","possible_defects")
#
# global const NN_positive1 = load("NN_positive_1_defects_N1000_W7.jld2","NN")
# global const possible_positive1_defects = load("NN_positive_1_defects_N1000_W7.jld2","possible_defects")
# global const NN_negative1 = load("NN_negative_1_defects_N1000_W7.jld2","NN")
# global const possible_negative1_defects = load("NN_negative_1_defects_N1000_W7.jld2","possible_defects")
