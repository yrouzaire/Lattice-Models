include("lattices.jl")
include("models.jl")
include("init_visu.jl")

include("core_methods.jl")
include("defects_methods.jl")

include("auxiliary.jl")
include("measurements.jl")

using BenchmarkTools, JLD2
global const WINDOW = 7
global const W21 = 2WINDOW+1

using Flux
using Flux:onecold, Chain, Dense, softmax
struct Reshape # has to be defined before Augmentor.jl is loaded to avoid conflicts
    shape
end
Reshape(args...) = Reshape(args)
(r::Reshape)(x) = reshape(x, r.shape)
Flux.@functor Reshape ()
# global const DAE_positive12 = load("DAE_positive12___03_11_2022.jld2","NN")

# global const NN_positive12 = load("NN_positive_12_defects_N1000_W7.jld2","NN")
# global const possible_positive12_defects = load("NN_positive_12_defects_N1000_W7.jld2","possible_defects")
# global const NN_negative12 = load("NN_negative_12_defects_N1000_W7.jld2","NN")
# global const possible_negative12_defects = load("NN_negative_12_defects_N1000_W7.jld2","possible_defects")
#
# global const NN_positive1 = load("NN_positive_1_defects_N1000_W7.jld2","NN")
# global const possible_positive1_defects = load("NN_positive_1_defects_N1000_W7.jld2","possible_defects")
# global const NN_negative1 = load("NN_negative_1_defects_N1000_W7.jld2","NN")
# global const possible_negative1_defects = load("NN_negative_1_defects_N1000_W7.jld2","possible_defects")
