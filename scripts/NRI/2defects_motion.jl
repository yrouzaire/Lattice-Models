cd("D:/Documents/Research/projects/LatticeModels")
    using DrWatson ; @quickactivate "LatticeModels"
    include(srcdir("LatticeModels.jl"))
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));

## Defect Tracking => (x,y,µ)(t) averaged over R ?

## Movies
r0 = round(Int,L/4)
tmax = 50 ; every = 0.5 ; times = 0:every:tmax
params_init["type2defect"] = [mus[i],mus[j]]
model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
# plot_thetas(thetas,model,lattice,defects=false)
update!(thetas,model,lattice,tmax=tmax)
