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

## Visual verification
model = VisionXY(params_phys,params_num)
thetas = init_thetas(lattice,init="hightemp",q=1,r0=60,float_type=float_type,type=["source","divergent"])
plot_theta(thetas,model,lattice)
tmax = 100
z = @elapsed update!(thetas,model,lattice,tmax)
    plot_theta(thetas,model,lattice)
    prinz(z)
    #ok tout semble correspondre à mes attentes
