"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Goal : Test
include("parameters.jl");
Ts = collect(1.8:0.1:2.2)
times_log = logspace(1,tmax,32)

polar_order   = zeros(length(Ts),length(times_log))
nematic_order = zeros(length(Ts),length(times_log))
C   = zeros(length(Ts),Int(L/2),length(times_log))
xi  = zeros(length(Ts),length(times_log))
n   = zeros(length(Ts),length(times_log))
thetas_save = zeros(Float16,length(Ts),length(times_log),L,L)

z = @elapsed for i in each(Ts)
    params["T"] = Ts[i]
    model   = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas  = init_thetas(lattice,params=params)

    token = 1
    while model.t < tmax
        update!(thetas,model,lattice)
        if model.t ≥ times_log[token]
            polar_order[i,token],nematic_order[i,token] = OP(thetas)
            correlation  = corr(thetas,model,lattice)
            C[i,:,token] = correlation
            xi[i,token]  = corr_length(correlation)
            n[i,token]   = number_defects(thetas,model,lattice)
            thetas_save[i,token,:,:] = thetas

            token = min(token+1,length(times_log))
        end
    end
end
prinz(z)

comments = "Goal: Recover TKT to check whether everything goes well. Model XY, on Triangular Lattice"
@save "data/TKT_$(symmetry)XY_r$(real).jld2" thetas_save polar_order nematic_order C xi n times_log Ts params runtime=z comments
