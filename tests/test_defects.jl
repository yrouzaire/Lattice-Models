include("../src/core_methods.jl");
include("../src/init_visu.jl");
include("../src/misc.jl");
using BenchmarkTools

## Parameters
include(srcdir("../parameters.jl"));

model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
# thetas = init_thetas(lattice,init="pair",q=1,r0=60,type=["source","divergent"])
spot_defects(thetas,model,lattice)
plot_theta(thetas,model,lattice,defects=true)

dft = DefectTracker(thetas,model,lattice)
total_tracked(dft)
defects_active(dft)
