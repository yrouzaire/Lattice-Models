using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


using BenchmarkTools
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)

## Tests get neighbours
# Physical Parameters
include(srcdir("../parameters.jl"));

## Benchmark update
model = MCXY(params)
lattice = TriangularLattice(L,periodic=true)
thetas = init_thetas(lattice,params=params)
update!(thetas,model,lattice)
@btime update!(thetas,model,lattice)
#= Runtimes
WN = Wrapped Normal ; VM = VonMises ; AF = antiferro
Avec WN, AF = false update takes 11 ms
Avec VM, AF = false update takes 23 ms
Avec WN, AF = true update takes 11 ms
Avec VM, AF = true update takes 25 ms
=#

## Tests update!()
model = XY(params_phys,params_num)
    @btime update!(thetas,model,lattice)

model = ForcedXY(params_phys,params_num)
    @btime update!(thetas,model,lattice)

model = VisionXY(params_phys,params_num)
    @btime update!(thetas,model,lattice)

model = MovingXY(params_phys,params_num)
thetas = init_thetas(model,lattice,init="2pair",q=1,r0=60,type=["source","divergent"])
update!(thetas,model,lattice)
@btime update!(thetas,model,lattice)
plot_theta(thetas,model,lattice)

# Tests update!(....,tmax)
tmax = 1
model = VisionXY(params_phys,params_num)
model.t
update!(thetas,model,lattice)
model.t
update!(thetas,model,lattice,tmax)
model.t

## Visual verification
include(srcdir("../parameters.jl"));
model = MCXY(params)
lattice = TriangularLattice(L,periodic=true)
thetas = init_thetas(lattice,params=params)
plot_thetas(thetas,model,lattice)
tmax = 10000
    z = @elapsed update!(thetas,model,lattice,tmax)
    plot_thetas(thetas,model,lattice)
    # ok tout semble correspondre à mes attentes
prinz(z)
plot_thetas(thetas,model,lattice,defects=true)

## Check chebychev Metric
model = ForceXY(params_phys,params_num)
lattice = SquareLattice(L,periodic=true,metric="chebychev")
thetas = init_thetas(model,lattice,init="hightemp",q=1,r0=60,float_type=float_type,type=["source","divergent"])
plot_theta(thetas,model,lattice)
tmax = 100
z = @elapsed update!(thetas,model,lattice,tmax)
plot_theta(thetas,model,lattice)
include("../src/defects_methods.jl")

plot_theta(thetas,model,lattice,defects=true)

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
