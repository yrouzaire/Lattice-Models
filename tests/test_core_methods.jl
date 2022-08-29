include("../src/core_methods.jl");
include("../src/init_visu.jl");

## Tests get neighbours
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

# Checks whether the number of NN is plausible
lattice = TriangularLattice(L)
i = 1; j=1
model = XY(params_phys,params_num)
    get_neighbours(model,lattice,i,j)
model = AXY(params_phys,params_num)
    get_neighbours(model,lattice,i,j)
model1 = VisionXY(params_phys,params_num)
    get_neighbours(model1,lattice,i,j,true)










##
