using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings,JLD2
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## XY Analysis
@load datadir("dft_XY.jld2") runtimes dfts R model lattice
histogram(runtimes/60,bins=round(Int,R/2))
MSD(dfts[1],model,lattice)
# TODO
t_bounds(dfts[1])
hasfield(typeof(model),:dt) ? dummy_dt = model.dt : dummy_dt = 1
creation_loc(dfts[1].defectsP[1])
creation_loc(dfts[1].defectsN[1])
