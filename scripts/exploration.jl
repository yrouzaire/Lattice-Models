cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


## Tracking defects over time (first single, then pair, then hightemp)
include(srcdir("../parameters.jl"));
    model = XY(params)
    lattice = TriangularLattice(L,periodic=true,single=true)
    thetas = init_thetas(lattice,params=params)

p=plot_thetas(randomly_rotate(thetas),model,lattice,defects=false)
    display_quiver!(p,randomly_rotate(thetas),5)
dft = DefectTracker(thetas,model,lattice)
dft.defectsP[1].type[1]
dft.defectsN[1].type[1]
update_and_track!(thetas,model,lattice,dft,100,5)
dft.defectsP[1].type


zoom_quiver(thetas,model,lattice,last_loc(dft.defectsP[1])...,9)
zoom_quiver(thetas,model,lattice,last_loc(dft.defectsN[1])...,9)
zoomN = zoom(thetas,lattice,last_loc(dft.defectsN[1])...,9)[2]
ind_type = onecold(NN(vec(zoomN)))
possible_defects[ind_type]




tmax,every = 100,10
z = @elapsed update_and_track_plot!(thetas,model,lattice,dft,tmax,every)
prinz(z)
dft

## Simple simulation and plots from hightemp
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = TriangularLattice(L,periodic=true,single=true)
thetas = init_thetas(lattice,params=params)
times = [200,300,500,800,1000]
times = [5,10,20,30,50,80,100]
# @elapsed z = update!(thetas,model,lattice,100)
token = 1
while model.t < times[end]
    update!(thetas,model,lattice,100)
    if model.t ≥ times[token]
        token = min(token+1,length(times))
        display(plot_thetas(thetas,model,lattice,defects=true))
    end
end
@elapsed z = update!(thetas,model,lattice,410)
plot_thetas(thetas,model,lattice)
spot_defects(thetas,model,lattice)
plot_thetas(thetas,model,lattice,defects=true)
