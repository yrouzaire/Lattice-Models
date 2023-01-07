cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));
using Random ; Random.seed!(1)

params["vision"] = 0.
model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,tmax=40)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)
thetasagg = zeros(Float32,2 .*size(thetas))
    thetasagg[1:L,1:L] = thetasagg[1:L,1+L:end] = thetasagg[1+L:end,1:L] = thetasagg[1+L:end,1+L:end] = thetas
    latticeagg = TriangularLattice(2L)
update!(thetasagg,model,latticeagg,tmax=50)  # updates until time = t
    p = plot_thetas(thetasagg,model,latticeagg,defects=false)

dft01 = DefectTracker(thetasagg,model,latticeagg,find_type=true)
histogram(last_types(dft),bins=50)
    histogram!(last_types(dft0),bins=50)
    histogram!(last_types(dft01),bins=50)
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

## Stability of µ = 0 et µ = pi for XY model, (Langevin et MonteCarlo)
include(srcdir("../parameters.jl"));
params_init["init"] = "single"
z = @elapsed for mu in [0,pi/2,pi,3pi/2]
    println("µ = $mu")
    params_init["type1defect"] = mu
    model = LangevinXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice,defects=false)
    Nsteps = 100 ; every = 0.5
    animation = @animate for i in 1:Nsteps
        update!(thetas,model,lattice,every)
        zoom_quiver(thetas,model,lattice,40,40,10)
    end
    mp4(animation,"films/soft_vision/stability_defects/stab_Langevin_µ$(round(mu,digits=2)).mp4")
end
prinz(z)

z = @elapsed for mu in 0#[0,pi/2,pi,3pi/2]
    println("µ = $mu")
    params_init["type1defect"] = mu
    model = MonteCarloXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice,defects=false)
    Nsteps = 100 ; every = 50
    animation = @animate for i in 1:Nsteps
        update!(thetas,model,lattice,every)
        zoom_quiver(thetas,model,lattice,40,40,10)
    end
    mp4(animation,"films/soft_vision/stability_defects/stab_MonteCarlo_µ$(round(mu,digits=2)).mp4")
end
prinz(z)

z = @elapsed for mu in [0,pi/2,pi,3pi/2]
    for vision in [0.2]
    println("µ = $(mu), vision = $vision")
    params_init["type1defect"] = mu
    params["vision"] = vision
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice,defects=false)
    Nsteps = 100 ; every = 0.5
    animation = @animate for i in 1:Nsteps
        update!(thetas,model,lattice,every)
        zoom_quiver(thetas,model,lattice,40,40,10)
    end
    mp4(animation,"films/soft_vision/stability_defects/stab_qneg_SoftVision_vision$(vision)_µ$(round(mu,digits=2)).mp4")
    end
end
prinz(z)

## µ(t) from µ(0) = 0 for SoftVisionXY
R = 2
visions = 1 .+ [0.02,0.05,0.1,0.2]
tmax = 10
for init_mu in [0,pi]
    for vision in visions
    for r in 1:R
        model = SoftVisionXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        dft = DefectTracker(thetas,model,lattice,find_type=true)
        update_and_track!(thetas,model,lattice,dft,tmax,every,find_type=true)
