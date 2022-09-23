"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Goal : Test
include("parameters.jl");
Ts = [0.1]
times_log = logspace(1,tmax,32)
times_lin = collect(transients:every:tmax)

polar_order   = zeros(length(Ts),length(times_log))
nematic_order = zeros(length(Ts),length(times_log))
C   = zeros(length(Ts),Int(L/2),length(times_log))
xi  = zeros(length(Ts),length(times_log))
n   = zeros(length(Ts),length(times_log))
thetas_save = zeros(Float16,length(Ts),length(times_log),L,L)

dfts = Vector{DefectTracker}(undef,length(Ts))

z = @elapsed for i in each(Ts)
    params["T"] = Ts[i]
    model   = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas  = init_thetas(model,lattice,params_init=params_init)

    token_log = 1 ; token_lin = 1
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ times_log[token_log]
            polar_order[i,token_log],nematic_order[i,token_log] = OP(thetas)
            correlation  = corr(thetas,model,lattice)
            C[i,:,token_log] = correlation
            xi[i,token_log]  = corr_length(correlation)
            n[i,token_log]   = number_defects(thetas,model,lattice)
            thetas_save[i,token_log,:,:] = thetas

            token_log = min(token_log+1,length(times_log))
        end
        if model.t ≥ times_lin[token_lin]
            if !isdefined(dfts,i) dfts[i] = DefectTracker(thetas,model,lattice,find_type=false)
            else update_DefectTracker!(dfts[i],thetas,model,lattice)
            end
            token_lin = min(token_lin+1,length(times_lin))
        end
    end
end
prinz(z)

comments = "Goal: Recover TKT to check whether everything goes well. Model XY, on Triangular Lattice"
@save "data/TKT_$(symmetry)XY_r$(real).jld2" dfts thetas_save polar_order nematic_order C xi n times_log times_lin Ts params runtime=z comments
