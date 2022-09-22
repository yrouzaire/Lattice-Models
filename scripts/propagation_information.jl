cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


## Code core function such as perturbation() or

## Test update 1D lattice
include(srcdir("../parameters.jl"));
lattice = Chain1D(L)
model = PropagationForcedXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice)

update!(thetas,model,lattice,300)
    plot_thetas(thetas,model,lattice)

## Test update 2D lattice
include(srcdir("../parameters.jl"));
lattice = SquareLattice(L)
model = PropagationForcedXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice)

update!(thetas,model,lattice,30)
    plot_thetas(thetas,model,lattice)
