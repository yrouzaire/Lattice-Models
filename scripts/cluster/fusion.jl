using JLD2
include("defects_methods.jl")

R = 40
Nb_thetas_save = 10
base_filename = "data/polarMovXY_rho1_A0"
indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
println("There are $(length(indices))/$R files.")

@load base_filename*"_r$(indices[1]).jld2" times_log times_lin Ts params comments

ns  = NaN*zeros(length(Ts),length(times_log),R)
Cs  = NaN*zeros(length(Ts),Int(params["L"]/2),length(times_log),R)
xis = NaN*zeros(length(Ts),length(times_log),R)
polar_orders = NaN*zeros(length(Ts),length(times_log),R)
nematic_orders = NaN*zeros(length(Ts),length(times_log),R)

dftss = Array{DefectTracker,2}(undef,length(Ts),R)

thetas_saves = NaN*zeros(Float16,length(Ts),length(times_log),params["L"],params["L"],Nb_thetas_save)
global token = 1
runtimes = NaN*zeros(R)

for r in indices
    @load base_filename*"_r$r.jld2" C n xi polar_order nematic_order dfts runtime
    ns[:,:,r] = n
    xis[:,:,r] = xi
    Cs[:,:,:,r] = C
    polar_orders[:,:,r] = polar_order
    nematic_orders[:,:,r] = nematic_order

    for i in 1:length(Ts)
        for n in 1:length(dfts[i].defectsP) dfts[i].defectsP[n].thetas_zoom = [zeros(Float32,2,2)] end
        for n in 1:length(dfts[i].defectsN) dfts[i].defectsN[n].thetas_zoom = [zeros(Float32,2,2)] end
    end
    dftss[:,r] = dfts

    if token ≤ Nb_thetas_save
        @load base_filename*"_r$r.jld2" thetas_save
        thetas_saves[:,:,:,:,token] = thetas_save
        global token += 1
    end
    runtimes[r] = runtime
end
@save base_filename*".jld2" times_log times_lin dftss Ts params comments polar_orders nematic_orders xis Cs ns thetas_saves runtimes R
println("Fusionned data saved in $(base_filename*".jld2") .")


# R = 10
# base_filename = "data/dft_XY"
# indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
# println("There are $(length(indices))/$R files.")
#
# # @unpack params,comments,model,lattice = JLD2.load(base_filename*"_r1.jld2")
# @load base_filename*"_r$(indices[1]).jld2" params comments model lattice
#
# dfts     = Vector{Any}(missing,R)
# runtimes = Vector{Any}(missing,R)
# for r in indices
#     @load base_filename*"_r$r.jld2" dft runtime
#     dfts[r] = dft
#     runtimes[r] = runtime
# end

# @save base_filename*".jld2" params comments model lattice dfts runtimes R
# println("Fusionned data saved in $(base_filename*".jld2") .")
