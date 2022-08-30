include("../src/core_methods.jl");
include("../src/init_visu.jl");
include("../src/misc.jl");

using BenchmarkTools

# Physical Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    Var = 0.1
    vision = π
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"vision"=>vision,"symmetry"=>symmetry)

# Numerical Parameters
dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

lattice = TriangularLattice(L,periodic=true,single=true)



get_vorticity(mod.(thetas,sym(model)),model,lattice,10,10)
