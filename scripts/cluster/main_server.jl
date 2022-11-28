"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

include("parameters_cluster.jl"); # imports Ts, As, rhos
times_log = logspace(1,tmax,10)
times_lin = collect(transients:every:tmax)

polar_order   = zeros(length(Ts),length(visions),length(rhos),length(times_log))
nematic_order = zeros(length(Ts),length(visions),length(rhos),length(times_log))
C   = zeros(length(Ts),length(visions),length(rhos),Int(L/2),length(times_log))
xi  = zeros(length(Ts),length(visions),length(rhos),length(times_log))
n   = zeros(length(Ts),length(visions),length(rhos),length(times_log))
thetas_save = zeros(Float16,length(Ts),length(visions),length(rhos),length(times_log),L,L)
# dfts = Array{Union{Missing,DefectTracker},3}(missing,length(Ts),length(As),length(rhos))

z = @elapsed for i in each(Ts) , j in each(As) , k in each(rhos)
    params["T"] = Ts[i]
    params["vision"] = visions[j]
    params["rho"] = rhos[k]

    model   = SoftVisionXY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas  = init_thetas(model,lattice,params_init=params_init)

    token_log = 1 ; token_lin = 1
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ times_log[token_log]
            polar_order[i,j,k,token_log],nematic_order[i,j,k,token_log] = OP(thetas)
            correlation  = corr(thetas,model,lattice)
            C[i,j,k,:,token_log] = correlation
            xi[i,j,k,token_log]  = corr_length(correlation)
            n[i,j,k,token_log]   = number_defects(thetas,model,lattice)
            thetas_save[i,j,k,token_log,:,:] = thetas

            token_log = min(token_log+1,length(times_log))
        end
        # if model.t ≥ times_lin[token_lin]
        #     if ismissing(dfts[i,j,k]) dfts[i,j,k] = DefectTracker(thetas,model,lattice,find_type=false)
        #     else update_DefectTracker!(dfts[i,j,k],thetas,model,lattice)
        #     end
        #     token_lin = min(token_lin+1,length(times_lin))
        # end
    end
end
prinz(z)

comments = "First simulations for SoftVisionXY to see what going on. A few visions, a few Ts. Model $(symmetry)SoftVisionXY, on Triangular Lattice"
@save "data/first_simus_NRI_r$(real).jld2" scanned_params thetas_save polar_order nematic_order C xi n times_log times_lin params runtime=z comments
# @save "data/$(symmetry)MovXY_r$(real).jld2" dfts thetas_save polar_order nematic_order C xi n times_log times_lin params runtime=z comments
