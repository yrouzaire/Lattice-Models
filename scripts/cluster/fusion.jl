using JLD2,Parameters

R = 40
base_filename = "data/dft_XY"
indices = [] ; for r in 1:R  if isfile(filename_base*"_r$r.jld") push!(indices,r) end end
println("There are $(length(indices))/$R files.")
@unpack params,comments,model,lattice = load(base_filename*"_r1.jld")

dfts     = Vector{Any}(missing,R)
runtimes = Vector{Any}(missing,R)
for r in indices
    dfts[r],runtimes[r] = load(base_filename*"_r$r.jld","dft","runtime")
end

@save base_filename*".jld" params comments model lattice dfts runtimes R
println("Fusionned data saved in $filename_fusion.")
