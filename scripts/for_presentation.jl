cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
fp = "plots/for_presentation_ForcedXY/"

#= In addition to the figures of the Frontiers in Physics, I need:
    * Photos pour les triangles, XY and Forced
    * Films avec pair de défaut
    *
    *
=#
include(srcdir("../parameters.jl"));
lattice = SquareLattice(L)
model   = XY(params)
thetas = thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice)
update!(thetas,model,lattice,tmax=100)
    plot_thetas(thetas,model,lattice,defects=false,colorbar=false)
savefig(fp*"photoFXY.png")



saving_times = 0:5:50
z = @elapsed anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients=0)
mp4(anim,)
