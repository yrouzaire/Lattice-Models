include("../src/core_methods.jl");
include("../src/init_visu.jl");
include("../src/misc.jl");
using BenchmarkTools

## Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    Var = 0.1
    vision = π
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"vision"=>vision,"symmetry"=>symmetry)

dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,init="isolated",q=1,type="source")
# thetas = init_thetas(lattice,init="pair",q=1,r0=60,type=["source","divergent"])
spot_defects(thetas,model,lattice)
plot_theta(thetas,model,lattice,defects=true)
