include("../src/core_methods.jl");
include("../src/init_visu.jl");

## Tests get neighbours
# Physical Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    params_phys = Dict("L"=>L,"T"=>T,"symmetry"=>symmetry)

# Numerical Parameters
dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

model = XY(params_phys,params_num)
lattice = TriangularLattice(L)
init_thetas!(model,lattice,init="hightemp")
plot_theta(model,lattice)

i = 1; j=1
get_neighbours(model,lattice,i,j)


L = 200
model = VisionXY(L,0.1,pi-0.1,"polar")
lattice = TriangularLattice(L)
init_thetas!(model,lattice,init="hightemp")
plot_theta(model,lattice)

i = 1; j=1
get_neighbours(model,lattice,i,j)

##
