cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Movies
include(srcdir("../parameters.jl"));

model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
saving_times = 0:100:20000 ; transients = saving_times[end]/2
z = @elapsed anim_extensile = movies(thetas,model,lattice,defects=false,saving_times=saving_times,transients=transients)
prinz(z)
mp4(anim_extensile,datadir("../films/test.mp4"))
mp4(anim_extensile,datadir("../films/active_extensile_A2.mp4"))
