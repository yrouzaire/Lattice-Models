"Copyright (c) 2022 Y.Rouzaire All Rights Reserved."

include("LatticeModels.jl") ;
using Plots # for plotting methods such as @animate etc
include("IDrealisation.jl") ;

## Scan parameters
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
z = @elapsed update!(thetas,model,lattice,10)
prinz(z)

base_filename = "data/test_save"
JLD.save(filename,"params",params,"runtime",z,"comments",comments,"thetas",thetas,"model",model,"lattice",lattice)
