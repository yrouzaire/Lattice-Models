using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Parameters
include(srcdir("../parameters.jl"));

## Benchmark update
model = XY(params_phys,params_num)
lattice = TriangularLattice(L,periodic=true)
thetas = init_thetas(model,lattice,init="hightemp",q=1,r0=60,float_type=float_type,type=["source","divergent"])
update!(thetas,model,lattice,1)
plot_theta(thetas,model,lattice)
