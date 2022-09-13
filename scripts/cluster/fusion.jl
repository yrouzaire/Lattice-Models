using JLD2
include("defects_methods.jl")

R = 40
base_filename = "data/dft_XY"
indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
println("There are $(length(indices))/$R files.")

# @unpack params,comments,model,lattice = JLD2.load(base_filename*"_r1.jld2")
@load base_filename*"_r$(indices[1]).jld2" params comments model lattice

dfts     = Vector{Any}(missing,R)
runtimes = Vector{Any}(missing,R)
for r in indices
    @load base_filename*"_r$r.jld2" dft runtime
    dfts[r] = dft
    runtimes[r] = runtime
end

@save base_filename*".jld2" params comments model lattice dfts runtimes R
println("Fusionned data saved in $(base_filename*".jld2") .")
