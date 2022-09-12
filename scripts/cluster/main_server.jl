"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots,JLD2 # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Scan parameters
include("parameters.jl");
model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
z = @elapsed update!(thetas,model,lattice,10)
prinz(z)

@save "data/test_save.jld2" params runtime=z comments thetas model lattice
