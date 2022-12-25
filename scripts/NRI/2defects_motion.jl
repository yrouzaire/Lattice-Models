cd("D:/Documents/Research/projects/LatticeModels")
    using DrWatson ; @quickactivate "LatticeModels"
    include(srcdir("LatticeModels.jl"))
    using Plots,ColorSchemes,LaTeXStrings
    gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));

## Visualise pair of defect
params_init["type2defect"] = [pi,pi]
    params_init["r0"] = 16
    params_init["phi"] = 0pi/2
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
    # plot_thetas(thetas,model,lattice,defects=true,title=titre)
    zoom_quiver(thetas,model,lattice,50,50,12)
    # zoom_quiver(thetas,model,lattice,50+round(Int,params_init["r0"]/2),50)
    title!(titre)
mod.(sum(params_init["type2defect"]) - params_init["phi"]+pi,2pi)
mod.(sum(params_init["type2defect"]) + params_init["phi"],2pi)


## Understand what is the µ of two defects
include(srcdir("../parameters.jl"));
    mus = 0:pi/8:2pi
    actual_mus = Vector{Tuple{Float64,Float64}}(undef,length(mus))
    for i in each(mus)
    r0 = round(Int,L/8)
    params_init["type2defect"] = [0,mus[i]]
    params_init["r0"] = r0
    params_init["phi"] = pi
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    dft = DefectTracker(thetas,model,lattice,find_type=true)
    # zoom_quiver(thetas,model,lattice,50,50,12,size=(700,700))
    plot_thetas(thetas,model,lattice,defects=true)

    actual_mus[i] = last_type(dft.defectsP[1]),last_type(dft.defectsN[1])
end
    p=plot(xlims=ylims=(0,2pi))
    for i in each(actual_mus)
        scatter!(actual_mus[i])
    end
    p
    plot!(x->mod(x-pi+2params_init["phi"],2pi))

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
    for p in 1#:1:length(phis)
        scatter!(actual_mus[p,:])
    end
    p
# find the x coordinate when y ~ 0
fphi = NaN*zeros(length(phis))
for p in each(phis)
    fphi[p] = mus[findmin([actual_mus[p,i][1] for i in each(mus)])[2]]
end
plot(phis,fphi)
    plot!(x->mod(pi+x,2pi))



## Defect Tracking => (x,y,µ)(t) averaged over R ?

## Movies
include(srcdir("../parameters.jl"));

tmax = 200 ; every = 1 ; times = 0:every:tmax
    params_init["r0"] = round(Int,32)
    mu_plus = pi
    mu_minus = pi/2
    params_init["type2defect"] = [mu_plus,mu_minus]
    params_init["phi"] = 0
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
    prinz(z)
    mp4(anim,"films/NRI/trajectories_2defects/µ$(mu_plus)_µ$(mu_minus)_phi$(round(params_init["phi"],digits=1)).mp4")

## Attraction ? Repulsion ? between a pair of defects
#=  Important control parameters of the dynamics of two defects:
    1. µ+ + µ- (the effective µs)
    2. phi (Free from  0 to 2pi)
To obtain a given µ+ + µ-, one has to set µ1 = 0 and µ2 = (µ+ + µ- + pi)/2
=#
include(srcdir("../parameters.jl"));
dµ = pi/8 ; sum_mus = Float32.(round.(collect(0:dµ:2pi),digits=3))
phis = sum_mus
r0 = round(Int,L/3)
tmax = 20
every = 1
R = 5
sigmas = [0.3]
rs = zeros(length(sum_mus),length(phis),length(sigmas),R)
dfts = Array{DefectTracker}(undef,length(sum_mus),length(phis),length(sigmas),R)
z = @elapsed for i in each(sum_mus) , j in each(phis)
    println((i-1)*length(sum_mus)+j,"/",length(sum_mus)*length(phis))
    for sig in each(sigmas)
    Threads.@threads for r in 1:R
        params["vision"] = sigmas[sig]
        params_init["type2defect"] = [0,(sum_mus[i]+pi)/2]
        params_init["phi"] = phis[j]
        params_init["r0"]  = r0
        model = SoftVisionXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        # plot_thetas(thetas,model,lattice,defects=false)
        update!(thetas,model,lattice,tmax=tmax)
        dft = DefectTracker(thetas,model,lattice)
        # update_and_track!(thetas,model,lattice,dft,tmax,every,find_type=true)
        # plot_thetas(thetas,model,lattice,defects=false)
        if number_active_defects(dft)>0
            rs[i,j,sig,r] = dist(lattice,last_loc(dft.defectsP[1]),last_loc(dft.defectsN[1]))
        else
            rs[i,j,sig,r] = 0
        end
        # dfts[i,j,sig,r] = dft
    end
    end
end
prinz(z)
rs_avg = nanmean(rs,4)[:,:,:,1]/r0
plot(xlabel=L"µ_{+} + µ_{-}",ylabel=L"ϕ",size=(470,400))
    heatmap!(sum_mus,phis,rs_avg[:,:,1],c=cgrad([:blue,:deepskyblue2,:white,:red,:red2]),aspect_ratio=1,colorbartitle="R(t)/R0",clims=(minimum(0),1))
    # plot!(mus,5pi/2 .+ pi/2 .- mus)
# savefig("plots/NRI/attraction_µµ_sigma0.34.png")
# @save "data/attraction_µµ_sigma0.3_zoom.jld2" mus rs rs_avg sigmas R tmax dµ every r0
# @load "data/attraction_µµ_sigma0.3.jld2" mus rs rs_avg sigmas R tmax dµ every r0

## Analyse d'images
plot_thetas(thetas,model,lattice,defects=true)
zoom_quiver(thetas,model,lattice,50,50,20,size=(700,700))
dft = DefectTracker(thetas,model,lattice,find_type=true)
dft.defectsN[1].pos[end]
tp = round(last_type(dft.defectsP[1]),digits=2)
tm = round(last_type(dft.defectsN[1]),digits=2)
