using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

#=Idea : for nematic/polar systems with densities 1 and below, find proportion of each defect type.
At hightemp, each defect should be equaprobable. If we find something very different after some relaxation time,
it means that the dynamics select some over the others.
Do the whole thing with L = 500 and trelax = 100 ?
Note 1 : with L = 300, polar symmetry and T=0.1,
    Langevin dynamics take t = 200 to arrive to 74 defects
    MCXY dynamics with proposal 2\sqrt T take t = 1500 to arrive to 76 defects
    MCXY dynamics with proposal 0-2pi take t = 6000 to arrive to also 76 defects
At first impression, the dynamics does not seem to change the nature of the defects
=#

## Code to zoom on the defects of a DefectTracker
n=0
n += 1
    d = dft.defectsP[n]
    loc = last_loc(d)
    zoom_quiver(thetas,model,lattice,loc...,10)
    title!(last_type(d))


# update_and_track!(thetas,model,lattice,dft,2100,5)
# plot_thetas(thetas,model,lattice,defects=false)
zoom_quiver(thetas,model,lattice,150,180,15)

## Generate configurations
include(srcdir("../parameters.jl"));
lattice = TriangularLattice(L)

# Results for MonteCarlo XY
params["symmetry"] = "nematic" ; params["rho"] = 1
model = MCXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice,tmax=3000)
number_defects(thetas,model,lattice)
plot_thetas(thetas,model,lattice,defects=false)
dft = DefectTracker(thetas,model,lattice,find_type=true)
number_defects_types(dft)*2/sum(number_defects_types(dft))


## Moving XY
params["symmetry"] = "nematic" ; params["rho"] = 0.95 ; params["A"] = 1.5
model = MovingXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice,tmax=3000)
    plot_thetas(thetas,model,lattice,defects=false)
number_defects(thetas,model,lattice)
dft = DefectTracker(thetas,model,lattice,find_type=true)
number_defects_types(dft)*2/sum(number_defects_types(dft))
