cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

## Visualise pair of defect
include(srcdir("../parameters.jl"));
    params_init["r0"] = 16
    params_init["mu_plus"] = 0
    params_init["mu_minus"] = nothing
    params_init["phi"] = pi/2
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    dft = DefectTracker(thetas,model,lattice,find_type=true)
    # zoom_quiver(thetas,model,lattice,50,50,12,size=(700,700))
    if number_active_defects(dft) > 0
        tp = string(round(last_type(dft.defectsP[1]),digits=2))
        tm = string(round(last_type(dft.defectsN[1]),digits=2))
        titre = L"µ_{+}="*tp*L" ; µ_{-}="*tm
    else
        titre = ""
    end
    plot_thetas(thetas,model,lattice,defects=false,title=titre)
    zoom_quiver(thetas,model,lattice,30,30,12)
    # zoom_quiver(thetas,model,lattice,50+round(Int,params_init["r0"]/2),50)
    title!(titre)

mod.(sum(params_init["type2defect"]) - params_init["phi"]+pi,2pi)
mod.(sum(params_init["type2defect"]) + params_init["phi"],2pi)


## Understand what is the µ of two defects
include(srcdir("../parameters.jl"));
    mus = 0:pi/8:2pi
    actual_mus = Matrix{Tuple{Float64,Float64}}(undef,length(mus),length(mus))
    for i in each(mus), j in each(mus)
    r0 = round(Int,L/8)
    params_init["type2defect"] = [0,mus[j]]
    params_init["r0"] = r0
    params_init["phi"] = pi/3
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    dft = DefectTracker(thetas,model,lattice,find_type=true)
    # zoom_quiver(thetas,model,lattice,50,50,12,size=(700,700))
    plot_thetas(thetas,model,lattice,defects=true)

    actual_mus[i,j] = last_type(dft.defectsP[1]),last_type(dft.defectsN[1])
end
    p=plot()
    for i in each(actual_mus)
        scatter!(actual_mus[i])
    end
    p
    plot!(x->mod(x-params_init["phi"],2pi))

## Find f(ϕ) such that µ+ - µ- = f(ϕ)
phis = 0:pi/16:2pi
mus = 0:pi/32:2pi
actual_mus = Matrix{Tuple{Float64,Float64}}(undef,length(phis),length(mus))
for p in each(phis) , i in each(mus)
        r0 = round(Int,L/8)
        params_init["type2defect"] = [0,mus[i]]
        params_init["r0"] = r0
        params_init["phi"] = phis[p]
        model = SoftVisionXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        dft = DefectTracker(thetas,model,lattice,find_type=true)
        actual_mus[p,i] = last_type(dft.defectsP[1]),last_type(dft.defectsN[1])
end
p=plot()
    for p in 1:1:length(phis)
        scatter!(actual_mus[p,:])
    end
    p
# find the x coordinate when y ~ 0
fphi = NaN*zeros(length(phis))
for p in each(phis)
    x = actual_mus[p,:]
    ind = nothing
    maxx = Inf
    for i in each(mus)
        if x[i][2] < maxx
            maxx = x[i][2]
            ind = i
        end
    end
    fphi[p] = mus[ind]
end
plot(phis,fphi)
    plot!(x->mod(-x,2pi))
