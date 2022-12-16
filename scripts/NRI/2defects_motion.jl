cd("D:/Documents/Research/projects/LatticeModels")
    using DrWatson ; @quickactivate "LatticeModels"
    include(srcdir("LatticeModels.jl"))
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));


## Understand what is the µ of two defects
r0 = round(Int,L/8)
    tmax = 50 ; every = 0.5 ; times = 0:every:tmax
    params_init["type2defect"] = [pi/2,pi]
    params_init["r0"] = r0
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    dft = DefectTracker(thetas,model,lattice,find_type=true)
    zoom_quiver(thetas,model,lattice,50,50,12,size=(700,700))
    # plot_thetas(thetas,model,lattice,defects=true)

last_type(dft.defectsP[1]),last_type(dft.defectsN[1])

## Defect Tracking => (x,y,µ)(t) averaged over R ?

## Movies
r0 = round(Int,L/4)
tmax = 50 ; every = 0.5 ; times = 0:every:tmax
params_init["type2defect"] = [0,0]
model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
plot_thetas(thetas,model,lattice,defects=false)
dft = DefectTracker(thetas,model,lattice,find_type=true)

update!(thetas,model,lattice,tmax=tmax)
z = @elapsed anim = @animate for tt in each(times)
    plot_thetas(thetas,model,lattice,defects=false)
    update!(thetas,model,lattice)
