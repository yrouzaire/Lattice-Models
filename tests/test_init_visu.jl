cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Tests Movies
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5)
model = MovingXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,init=init,q=1,r0=60,type=["source","divergent"])
    plot_theta(thetas,model,lattice)
saving_times = 0:100:20000 ; transients = Inf
# update!(thetas,model,lattice,100)
# plot_theta(thetas,model,lattice,defects=true)
z = @elapsed anim = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
prinz(z)
mp4(anim,"D:/Documents/Research/projects/lattice_models/Lattice Models/films/test3.mp4")

## Basic Tests
thetas = init_thetas(model,lattice,init="isolated",q=1,type="source")
    plot_theta(thetas,model,lattice)
&

thetas = init_thetas(model,lattice,init="isolated",q=-1,type="convergent")
    plot_theta(thetas,model,lattice)
&

thetas = init_thetas(model,lattice,init="isolated",q=1/2,type="source")
    plot_theta(thetas,model,lattice)
&

thetas = init_thetas(model,lattice,init="pair",q=1,r0=60,type=["source","divergent"])
    plot_theta(thetas,model,lattice)
&

thetas = init_thetas(model,lattice,init="2pair",q=1,r0=60,type=["source","divergent"])
    plot_theta(thetas,model,lattice)
&

thetas = init_thetas(model,lattice,init="pair",q=1/2,r0=60,type=["source","divergent"])
    plot_theta(thetas,model,lattice)
&

thetas = init_thetas(model,lattice,init="2pair",q=1/2,r0=60,type=["source","divergent"])
    plot_theta(thetas,model,lattice)

## Test init pair the new way
include(srcdir("../parameters.jl"));
    params["init"] = "pair"
    params["type2defect"] = "pair1"

    model = XY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(lattice,params=params)

p = plot_thetas(thetas,model,lattice,colorbar=false)
update!(thetas,model,lattice,10)
window = 5
    display_quiver!(p,thetas,window)
    xlims!(1,2window+1) ; ylims!(1,2window+1)

## Test the zoom with periodic lattice
include(srcdir("../parameters.jl"));

model = XY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(lattice,params=params)
p = plot_thetas(thetas,model,lattice,colorbar=false)
update!(thetas,model,lattice,50)
p = plot_thetas(thetas,model,lattice,colorbar=false)

window = 30
thetas_zoom = zoom(thetas,lattice,98,96,window)[2]
    p = plot_thetas(thetas_zoom,model,lattice,colorbar=false)

window = 12
    zoom_quiver(thetas,model,lattice,98,96,window)
