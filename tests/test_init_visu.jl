include("../src/init_visu.jl");

# Physical Parameters
L = 300
    T = 0.1
    symmetry = "polar"
    Var = 0.1
    vision = π
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"vision"=>vision,"symmetry"=>symmetry)

# Numerical Parameters
dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

model = XY(params_phys,params_num)
lattice = TriangularLattice(L)
init_thetas!(model,lattice,init="2pair",q=1,r0=60,type=["source","divergent"])
    plot_theta(model,lattice)
&

## Basic Tests
init_thetas!(model,lattice,init="isolated",q=1,type="source")
    plot_theta(model,lattice)
&

init_thetas!(model,lattice,init="isolated",q=-1,type="convergent")
    plot_theta(model,lattice)
&

init_thetas!(model,lattice,init="isolated",q=1/2,type="source")
    plot_theta(model,lattice)
&

init_thetas!(model,lattice,init="pair",q=1,r0=60,type=["source","divergent"])
    plot_theta(model,lattice)
&

init_thetas!(model,lattice,init="2pair",q=1,r0=60,type=["source","divergent"])
    plot_theta(model,lattice)
&

init_thetas!(model,lattice,init="pair",q=1/2,r0=60,type=["source","divergent"])
    plot_theta(model,lattice)
&

init_thetas!(model,lattice,init="2pair",q=1/2,r0=60,type=["source","divergent"])
    plot_theta(model,lattice)
&

## Tests vortices TODO
