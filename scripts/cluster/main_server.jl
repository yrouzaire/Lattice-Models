"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Goal : compare motion of a defect pair
include("parameters.jl");
model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
dft = DefectTracker(thetas,model,lattice)

tmax,every = 1000,5
z = @elapsed update_and_track!(thetas,model,lattice,dft,tmax,every)
prinz(z)

comments = ""

@save "data/dft_XY_r$(real).jld2" params runtime=z comments thetas model lattice dft tmax every
