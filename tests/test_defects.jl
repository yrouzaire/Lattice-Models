include("../src/core_methods.jl");
include("../src/init_visu.jl");
include("../src/misc.jl");
using BenchmarkTools

## Parameters
include(srcdir("../parameters.jl"));

lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,init="isolated",q=1,type="source")
# thetas = init_thetas(lattice,init="pair",q=1,r0=60,type=["source","divergent"])
spot_defects(thetas,model,lattice)
plot_theta(thetas,model,lattice,defects=true)
