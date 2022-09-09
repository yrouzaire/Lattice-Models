using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl")) ;
using BenchmarkTools,Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Code for update_and_track
include(srcdir("../parameters.jl"));
    model = MovingXY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    dft = DefectTracker(thetas,model,lattice)

    # update_and_track!(thetas,model,lattice,dft,50,.1)
    update_and_track_plot!(thetas,model,lattice,dft,500,10,defects=true)

number_active_defects(dft)
last_loc(dft.defectsN[1])
last_loc(dft.defectsP[1])
dft
plot_thetas(thetas,model,lattice,defects=true)
&

## Preconditionning of the theta field
include(srcdir("../parameters.jl"));
    model = MovingXY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    # update!(thetas,model,lattice,100)
    plot_thetas(thetas,model,lattice,defects=false)
    # thetas0 = copy(thetas)
    # plot_thetas(precondition!(thetas0,model,lattice),model,lattice,defects=true)

















&
