include("../src/core_methods.jl");
include("../src/init_visu.jl");
include("../src/misc.jl");

using BenchmarkTools
## Tests get neighbours
# Physical Parameters
L = 200
    T = 0.1
    symmetry = "polar"
    Var = 0.1
    A = 0.5
    vision = 4π/3
    rho = 1
    params_phys = Dict("L"=>L,"T"=>T,"Var"=>Var,"A"=>A,"rho"=>rho,"vision"=>vision,"symmetry"=>symmetry)

# Numerical Parameters
dt = 1E-2
    float_type = Float32
    params_num  = Dict("dt"=>dt,"float_type"=>float_type)

lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(model,lattice,init="2pair",q=1,r0=60,type=["source","divergent"])

## Checks whether the number of NN is plausible
i = 10; j=10
model = XY(params_phys,params_num)
    @btime get_neighbours(thetas,model,lattice,i,j)

model = ForcedXY(params_phys,params_num)
    @btime get_neighbours(thetas,model,lattice,i,j)

model = MovingXY(params_phys,params_num)
    @btime get_neighbours(thetas,model,lattice,i,j)

model = VisionXY(params_phys,params_num)
angle_neighbours = get_neighbours(thetas,model,lattice,i,j,true)
sum_influence_neighbours(thetas[i,j],angle_neighbours,model,lattice)
@btime get_neighbours(thetas,model,lattice,i,j,true)
@btime sum_influence_neighbours(thetas[i,j],angle_neighbours,model,lattice)


## Tests update!()
model = XY(params_phys,params_num)
    @btime update!(thetas,model,lattice)
model = ForcedXY(params_phys,params_num)
    @btime update!(thetas,model,lattice)
model = VisionXY(params_phys,params_num)
    @btime update!(thetas,model,lattice)

# Tests update!(....,tmax)
tmax = 1
model = VisionXY(params_phys,params_num)
model.t
update!(thetas,model,lattice)
model.t
update!(thetas,model,lattice,tmax)
model.t

## Visual verification
model = VisionXY(params_phys,params_num)
lattice = TriangularLattice(L,periodic=true)
lattice = SquareLattice(L,periodic=true)
thetas = init_thetas(lattice,init="hightemp",q=1,r0=60,float_type=float_type,type=["source","divergent"])
plot_theta(thetas,model,lattice)
tmax = 100
z = @elapsed update!(thetas,model,lattice,tmax)
    plot_theta(thetas,model,lattice)
    # ok tout semble correspondre à mes attentes
prinz(z)

## Check chebychev Metric
model = ForceXY(params_phys,params_num)
lattice = SquareLattice(L,periodic=true,metric="chebychev")
thetas = init_thetas(lattice,init="hightemp",q=1,r0=60,float_type=float_type,type=["source","divergent"])
plot_theta(thetas,model,lattice)
tmax = 100
z = @elapsed update!(thetas,model,lattice,tmax)
plot_theta(thetas,model,lattice)
include("../src/defects_methods.jl")

plot_theta(thetas,model,lattice,defects=true)
