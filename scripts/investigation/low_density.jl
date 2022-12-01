cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));
model = SPP(params) ; lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
tmax = 1000
update!(thetas,model,lattice,tmax=tmax)
    plot_thetas(thetas,model,lattice,defects=false)
delta = 100
update!(thetas,model,lattice,delta)
    plot_thetas(thetas,model,lattice,defects=false)

zoom_quiver(thetas,model,lattice,165,85)
