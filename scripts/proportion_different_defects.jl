using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

#=Idea : for nematic/polar systems with densities 1 and below, find proportion of each defect type.
At hightemp, each defect should be equaprobable. If we find something very different after some relaxation time,
it means that the dynamics select some over the others.
Do the whole thing with L = 500 and trelax = 100 ? =#
include(srcdir("../parameters.jl"));
lattice = TriangularLattice(L)

# Results for polar symmetry and rho = 1
params["symmetry"] = "polar" ; params["rho"] = 1
model = XY(params)
thetas = init_thetas(lattice,params=params)
@time update!(thetas,model,lattice,100)
plot_thetas(thetas,model,lattice,defects=false)
dft = DefectTracker(thetas,model,lattice,find_type=true)

number_defects_types(dft)*2/sum(number_defects_types(dft))

# Results for nematic symmetry and rho = 1
params["symmetry"] = "nematic" ; params["rho"] = 1
model = XY(params)
thetas = init_thetas(lattice,params=params)
@time update!(thetas,model,lattice,100)
plot_thetas(thetas,model,lattice,defects=false)
dft = DefectTracker(thetas,model,lattice,find_type=true)
number_defects_types(dft)*2/sum(number_defects_types(dft))

## Moving XY
# Results for nematic symmetry and rho = 1 and A = 0
params["symmetry"] = "nematic" ; params["rho"] = 0.9 ; params["A"] = 0
model = MovingXY(params)
thetas = init_thetas(lattice,params=params)
@time update!(thetas,model,lattice,1000)
plot_thetas(thetas,model,lattice,defects=false)
dft = DefectTracker(thetas,model,lattice,find_type=true)
number_defects_types(dft)*2/sum(number_defects_types(dft))

zoom_quiver(thetas,model,lattice,41,65,9)


# Results for nematic symmetry and rho < 1
