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
global const DAE_positive1  = load("DAE_positive1___15_12_2022.bson","DAE")
global const DAE_negative1  = load("DAE_negative1___15_12_2022.bson","DAE")
# global const DAE_positive12 = load("DAE_positive12___15_12_2022.bson","DAE")
# global const DAE_negative12 = load("DAE_negative12___15_12_2022.bson","DAE")
