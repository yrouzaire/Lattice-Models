using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));

#=Idea : for nematic/polar systems with densities 1 and below, find proportion of each defect type.
At hightemp, each defect should be equaprobable. If we find something very different after some relaxation time,
it means that the dynamics select some over the others.
Do the whole thing with L = 500 and trelax = 100 ? =#

# Results for polar symmetry and rho = 1
# Results for polar symmetry and rho < 1
# Results for nematic symmetry and rho = 1
# Results for nematic symmetry and rho < 1
