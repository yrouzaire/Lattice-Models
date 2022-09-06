using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl")) ;
using BenchmarkTools,Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Code for update_and_track
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)

dft = DefectTracker(thetas,model,lattice)
total_tracked(dft)
defects_active(dft)
