using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl")) ;
using BenchmarkTools,Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Code for update_and_track
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = SquareLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
dft = DefectTracker(thetas,model,lattice)

update_and_track!(thetas,model,lattice,dft,5,.1)
update_and_track_plot!(thetas,model,lattice,dft,10,1,defects=true)


MSD(dft,model,lattice)
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

## Issue locating defects v2
include(srcdir("../parameters.jl"));
    params["symmetry"] = "nematic"
    params["init"] = "single"
    params["q"] = 1/2
    params["rho"] = 1
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    # update!(thetas,model,lattice,1)
    p = plot_thetas(thetas,model,lattice,defects=true)


## Issue locating defects v1
include(srcdir("../parameters.jl"));
    params["symmetry"] = "polar"
    params["q"] = 1
    model = XY(params)
    lattice = SquareLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)
    z = @elapsed update!(thetas,model,lattice,0.3)

thetas_zoom = thetas[15:35,40:60]
    p=plot_thetas(thetas_zoom,model,lattice,defects=true)
    window = 10
    display_quiver!(p,thetas_zoom,window)
    xlims!(1,2window+1) ; ylims!(1,2window+1)
&




&
cols = cgrad([:black,:blue,:green,:orange,:red,:black])
heatmap(mod.(thetas,2pi)')
heatmap(thetas')
arclength(-3.1415,3.1415,pi)
arclength(3.1415,.0,pi)
plot_thetas(thetas,model,lattice,defects=true)
spot_defects(thetas,model,lattice)










&
