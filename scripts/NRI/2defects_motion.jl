cd("D:/Documents/Research/projects/LatticeModels")
    using DrWatson ; @quickactivate "LatticeModels"
    include(srcdir("LatticeModels.jl"))
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));


## Understand what is the µ of two defects
include(srcdir("../parameters.jl"));
mus = 0:pi/8:2pi
actual_mus = Matrix{Tuple{Float64,Float64}}(undef,length(mus),length(mus))
for i in each(mus), j in each(mus)
    r0 = round(Int,L/8)
    tmax = 10 ; every = 0.5 ; times = 0:every:tmax
    params_init["type2defect"] = [mus[i],mus[j]]
    params_init["r0"] = r0
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    dft = DefectTracker(thetas,model,lattice,find_type=true)
    # zoom_quiver(thetas,model,lattice,50,50,12,size=(700,700))
    # plot_thetas(thetas,model,lattice,defects=true)

    actual_mus[i,j] = last_type(dft.defectsP[1]),last_type(dft.defectsN[1])
end
p=plot()
    for i in each(actual_mus)
        scatter!(actual_mus[i])
    end
    plot!(x->mod(x-2,2pi))
    p
## Defect Tracking => (x,y,µ)(t) averaged over R ?

## Movies
tmax = 50 ; every = 0.25 ; times = 0:every:tmax
    params_init["r0"] = round(Int,32)
    mu_plus = 0
    mu_minus = 0
    params_init["type2defect"] = [mu_plus,mu_minus]
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice,defects=false)
    # dft = DefectTracker(thetas,model,lattice,find_type=true)

    z = @elapsed anim = @animate for tt in each(times)
    dft = DefectTracker(thetas,model,lattice,find_type=true)
    if number_active_defects(dft) > 0
        tp = string(round(last_type(dft.defectsP[1]),digits=2))
        tm = string(round(last_type(dft.defectsN[1]),digits=2))
        titre = L"µ_{+}="*tp*L" ; µ_{-}="*tm
    else
        titre = ""
    end

    p=plot_thetas(thetas,model,lattice,defects=true,size=(512,512),title=titre)
    update!(thetas,model,lattice,every)
    p
end
    mp4(anim,"films/NRI/trajectories_2defects/4defects_µ$(mu_plus)_µ$(mu_minus).mp4")
prinz(z)

## Analyse d'images
plot_thetas(thetas,model,lattice,defects=true)
zoom_quiver(thetas,model,lattice,50,50,20,size=(700,700))
dft = DefectTracker(thetas,model,lattice,find_type=true)
dft.defectsN[1].pos[end]
tp = round(last_type(dft.defectsP[1]),digits=2)
tm = round(last_type(dft.defectsN[1]),digits=2)
