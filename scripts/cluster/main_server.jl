"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Attraction ? Repulsion ? between a pair of defects
include("parameters_cluster.jl");
dµ = pi/2 ; mus = Float32.(round.(collect(0:dµ:2π),digits=2))
r0 = round(Int,L/3)
tmax = 100 ; every = 5 ; times = every:every:tmax
sigmas = [0.1,0.3]
separations = zeros(length(mus),length(mus),length(sigmas))
dfts = Array{DefectTracker}(undef,length(mus),length(mus),length(sigmas))
z = @elapsed for i in each(mus) , j in each(mus)
    # println((i-1)*length(mus)+j,"/",length(mus)^2)
    for sig in each(sigmas)
        params["vision"] = sigmas[sig]
        params_init["type2defect"] = [mus[i],mus[j]]
        model = SoftVisionXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        # plot_thetas(thetas,model,lattice,defects=false)
        token = 1
        while model.t < tmax
            update!(thetas,model,lattice)
            if model.t ≥ times[token]
                dft = DefectTracker(thetas,model,lattice)
                if number_active_defects(dft) > 0
                    separations[i,j,sig,token] = dist(lattice,last_loc(dft.defectsP[1]),last_loc(dft.defectsN[1]))
                else
                    separations[i,j,sig,token] = 0
                end
                token = min(token+1,length(times))
            end
        end
    end
end
prinz(z)
comments = "Manually create a config with two defects µ+ and µ- and let them evolve in time. Monitor their separation R(t) as a function of µ+ µ-"
@save "data/separations_µµ_r$(real).jld2" mus dµ separations r0 sigmas tmax every times params runtime=z comments

## Systemic probes
# include("parameters_cluster.jl");
# polar_order   = zeros(length(Ts),length(visions),length(times_log))
# nematic_order = zeros(length(Ts),length(visions),length(times_log))
# C   = zeros(length(Ts),length(visions),Int(L/2),length(times_log))
# xi  = zeros(length(Ts),length(visions),length(times_log))
# n   = zeros(length(Ts),length(visions),length(times_log))
# thetas_save = zeros(Float16,length(Ts),length(visions),length(times_log),L,L)
# # dfts = Array{Union{Missing,DefectTracker},3}(missing,length(Ts),length(visions))
#
# z = @elapsed for i in each(Ts) , j in each(visions)
#     params["T"] = Ts[i]
#     params["vision"] = visions[j]
#
#     model   = SoftVisionXY(params)
#     lattice = TriangularLattice(L,periodic=true,single=true)
#     thetas  = init_thetas(model,lattice,params_init=params_init)
#
#     token_log = 1 ; token_lin = 1
#     while model.t < tmax
#         update!(thetas,model,lattice)
#         if model.t ≥ times_log[token_log]
#             polar_order[i,j,token_log],nematic_order[i,j,token_log] = OP(thetas)
#             correlation  = corr(thetas,model,lattice)
#             C[i,j,:,token_log] = correlation
#             xi[i,j,token_log]  = corr_length(correlation)
#             n[i,j,token_log]   = number_defects(thetas,model,lattice)
#             thetas_save[i,j,token_log,:,:] = thetas
#
#             token_log = min(token_log+1,length(times_log))
#         end
#         # if model.t ≥ times_lin[token_lin]
#         #     if ismissing(dfts[i,j,k]) dfts[i,j,k] = DefectTracker(thetas,model,lattice,find_type=false)
#         #     else update_DefectTracker!(dfts[i,j,k],thetas,model,lattice)
#         #     end
#         #     token_lin = min(token_lin+1,length(times_lin))
#         # end
#     end
# end
# prinz(z)
#
# comments = "First simulations for SoftVisionXY to see what going on. A few visions, a few Ts. Model $(symmetry)SoftVisionXY, on Triangular Lattice"
# @save "data/first_simus_NRI_r$(real).jld2" scanned_params thetas_save polar_order nematic_order C xi n times_log times_lin params runtime=z comments
# @save "data/$(symmetry)MovXY_r$(real).jld2" dfts thetas_save polar_order nematic_order C xi n times_log times_lin params runtime=z comments
