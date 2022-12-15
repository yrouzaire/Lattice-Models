using JLD2
include("defects_methods.jl")

# R = 10
# Nb_thetas_save = 10
# base_filename = "data/first_simus_NRI"
# indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
# println("There are $(length(indices))/$R files.")

# @load base_filename*"_r$(indices[1]).jld2" times_log times_lin params scanned_params comments
# Ts = scanned_params["Ts"] ; visions = scanned_params["visions"] ; L = params["L"]
#
# ns  = NaN*zeros(length(Ts),length(visions),length(times_log),R)
# Cs  = NaN*zeros(length(Ts),length(visions),Int(L/2),length(times_log),R)
# xis = NaN*zeros(length(Ts),length(visions),length(times_log),R)
# polar_orders = NaN*zeros(length(Ts),length(visions),length(times_log),R)
# nematic_orders = NaN*zeros(length(Ts),length(visions),length(times_log),R)
#
# # dftss = Array{Union{Missing,DefectTracker},4}(missing,length(Ts),length(visions),R)
#
# thetas_saves = NaN*zeros(Float16,length(Ts),length(visions),length(times_log),L,L,Nb_thetas_save)
# global token = 1
# runtimes = NaN*zeros(R)
#
# for r in indices
#     @load base_filename*"_r$r.jld2" C n xi polar_order nematic_order runtime
#     ns[:,:,:,r] = n
#     xis[:,:,:,r] = xi
#     Cs[:,:,:,:,r] = C
#     polar_orders[:,:,:,r] = polar_order
#     nematic_orders[:,:,:,r] = nematic_order
#
#     # for i in each(dfts)
#     #     for n in 1:length(dfts[i].defectsP) dfts[i].defectsP[n].thetas_zoom = [zeros(Float32,2,2)] end
#     #     for n in 1:length(dfts[i].defectsN) dfts[i].defectsN[n].thetas_zoom = [zeros(Float32,2,2)] end
#     # end
#     # dftss[:,r] = dfts
#
#     if token ≤ Nb_thetas_save
#         @load base_filename*"_r$r.jld2" thetas_save
#         thetas_saves[:,:,:,:,:,token] = thetas_save
#         global token += 1
#     end
#     runtimes[r] = runtime
# end
# @save base_filename*".jld2" scanned_params times_log times_lin params comments polar_orders nematic_orders xis Cs ns thetas_saves runtimes R
# # @save base_filename*".jld2" scanned_params times_log times_lin dfts=dftss params comments polar_orders nematic_orders xis Cs ns thetas_saves runtimes R
# println("Fusionned data saved in $(base_filename*".jld2") .")


R = 40
base_filename = "data/separations_µµ"
indices = [] ; for r in 1:R  if isfile(base_filename*"_r$r.jld2") push!(indices,r) end end
println("There are $(length(indices))/$R files.")

@load base_filename*"_r$(indices[1]).jld2" mus separations sigmas tmax every times dµ r0 params comments

separations = NaN*zeros(length(mus),length(mus),length(sigmas),R)
runtimes = Vector{Any}(missing,R)
for r in indices
    @load base_filename*"_r$r.jld2" runtime separations
    separationss[:,:,:,r] = separations
    runtimes[r] = runtime
end

@save base_filename*".jld2" mus separationss sigmas tmax every times dµ r0 params comments runtimes R
println("Fusionned data saved in $(base_filename*".jld2") .")
