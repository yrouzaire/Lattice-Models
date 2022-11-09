cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));
using Random ; Random.seed!(1)

params["vision"] = 0.2
model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,tmax=30)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)
thetasagg = zeros(Float32,2 .*size(thetas))
    thetasagg[1:L,1:L] = thetasagg[1:L,1+L:end] = thetasagg[1+L:end,1:L] = thetasagg[1+L:end,1+L:end] = thetas
    latticeagg = TriangularLattice(2L)
update!(thetasagg,model,latticeagg,tmax=50)  # updates until time = t
    p = plot_thetas(thetasagg,model,latticeagg,defects=false)

## Movies
include(srcdir("../parameters.jl"));
params["vision"] = 0.2
# params["symmetry"] = "polar"
# params_init["init"] = "single"
params_init["type2defect"] = "pair2"
model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice,defects=false)

every = 0.5 ; tmax = 40 ; transients = Inf#round(Int,tmax*0.8) # defects are not plotted before t ≥ transients (if defects=true)
saving_times = every:every:tmax
z = @elapsed animation = movies(thetas,model,lattice,defects=true,saving_times=saving_times,transients=transients)
prinz(z)
filename = "films/soft_vision/single_µpi_$(params["symmetry"])_vision_$(model.vision).mp4"
mp4(animation,filename,fps=20)

plot_thetas(thetas,model,lattice,defects=false)
zoom_quiver(thetas,model,lattice,123,22,12)
