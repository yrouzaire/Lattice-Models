"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Goal : scan
include("parameters.jl");
model   = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas  = init_thetas(lattice,params=params)
dft     = DefectTracker(thetas,model,lattice)

tmax,every = Int(500),10
z = @elapsed while model.t < tmax
    update!(thetas,model,lattice,model.t + every)
    update_DefectTracker!(dft,thetas,model,lattice)
    push!(OPs,OP(thetas))
    push!(,)
end
prinz(z)

comments = ""

@save "data/dft_XY_r$(real).jld2" params runtime=z comments thetas model lattice dft tmax every
