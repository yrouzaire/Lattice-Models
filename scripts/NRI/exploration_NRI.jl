cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));

## Maximal sigma (non reciprocity)
include(srcdir("../parameters.jl"));
params["vision"] = 0.1
model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    # lattice = SquareLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    update!(thetas,model,lattice,100)
    plot_thetas(thetas,model,lattice,defects=false)


## Attraction ? Repulsion ? between a pair of defects
include(srcdir("../parameters.jl"));
dµ = pi/32 ; mus = Float32.(round.(collect(0:dµ:2pi),digits=3))
r0 = round(Int,L/3)
tmax = 50
every = 1
R = 10
sigmas = [0.3]
rs = zeros(length(mus),length(mus),length(sigmas),R)
dfts = Array{DefectTracker}(undef,length(mus),length(mus),length(sigmas),R)
z = @elapsed for i in each(mus) , j in each(mus)
    println((i-1)*length(mus)+j,"/",length(mus)^2)
    for sig in each(sigmas)
    Threads.@threads for r in 1:R
        params["vision"] = sigmas[sig]
        params_init["type2defect"] = [mus[i],mus[j]]
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
rs_avg = nanmean(rs,5)[:,:,:,:,1]/r0
plot(xlabel=L"µ_{+}",ylabel=L"µ_{-}",size=(470,400))
    heatmap!(mus,mus,rs_avg[:,:,2,end],c=cgrad([:blue,:deepskyblue2,:white,:red,:red2]),aspect_ratio=1,colorbartitle="R(t)/R0",clims=(minimum(0),1))
    # plot!(mus,5pi/2 .+ pi/2 .- mus)
# savefig("plots/NRI/attraction_µµ_sigma0.34.png")
# @save "data/attraction_µµ_sigma0.3_zoom.jld2" mus rs rs_avg sigmas R tmax dµ every r0
# @load "data/attraction_µµ_sigma0.3.jld2" mus rs rs_avg sigmas R tmax dµ every r0

## Stability of µ over time
tmax = 10 ; every = 0.1 ; times = collect(0:every:tmax+every)
ϵ = 0.
µ0s = [0,pi/2,pi,3pi/2] .+ ϵ
reals = 50
visions = [0,0.05,0.1,0.2,0.4]
mus = zeros(length(µ0s),length(visions),length(times)-1,reals)
z = @elapsed for i in each(µ0s)
    for j in each(visions)
        Threads.@threads for r in 1:reals
            println("r = $r ; µ0 = ",round(µ0s[i],digits=1)," ; vision = ",visions[j])
            params["L"] = 200 ; params["T"] = 0.1
            params["symmetry"] = "polar" ; params["rho"] = 1 ; params["vision"] = visions[j]
            params_init["init"] = "single" ; params_init["type1defect"] =  µ0s[i] # initial µ
            params_init["q"] = 1
            model = SoftVisionXY(params)
            lattice = TriangularLattice(L)
            thetas = init_thetas(model,lattice,params_init=params_init)

            dft = DefectTracker(thetas,model,lattice,find_type=true)
            update_and_track!(thetas,model,lattice,dft,tmax,every,find_type=true)
            mus[i,j,:,r] = dft.defectsP[1].type
        end
    end
end
    prinz(z) # 5 minutes for a lot of curves, very fast
mus_avg = mean(mus,dims=4)[:,:,:,1]
p=plot(ylims=(0,2pi),xlabel="t",ylabel="µ(t)")
    for i in [3,4]#each(µ0s)
        for j in each(visions)
            # for r in 1:reals
            #     plot!(times,smooth(mus[i,j,:,r],over=1),m=true,line=false,c=i)
            # end
            plot!(times[1:end-1],mus_avg[i,j,:],c=j,lw=2)#,label="σ = $(visions[j])",rib=0)
        end
    end
    p
savefig("plots/µt.png")
## Movies
include(srcdir("../parameters.jl"));

lattice = TriangularLattice(L)
model = SoftVisionXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
saving_times = 0:0.5:200 ; transients = saving_times[end]/2
z = @elapsed anim = movies(thetas,model,lattice,defects=false,saving_times=saving_times,transients=transients)
prinz(z)
    mp4(anim,datadir("../films/soft_vision/NRItest.mp4"))

## Movies and save configs
using Random
my_seed = 2
Random.seed!(my_seed)
include(srcdir("../parameters.jl"));

lattice = TriangularLattice(L)
model = SoftVisionXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
saving_times = 0:1:300; transients = saving_times[end]/2
thetas_saves = zeros(Float16,L,L,length(saving_times))
token = 1
z = @elapsed while model.t < saving_times[end]
    update!(thetas,model,lattice)
    if model.t ≥ saving_times[token]
        println("t = $(round(model.t,digits=1))")
        thetas_saves[:,:,token] = thetas
        token = min(token+1,length(saving_times))
    end
end
@save "data/etude_films/NRI_vision$(vision)_T$(T)_seed$(my_seed).jld2" vision T L seed=my_seed thetas_saves model lattice saving_times
# movie of all system
anim = @animate for tt in 1:length(saving_times)
    plot_thetas(thetas_saves[:,:,tt],model,lattice,size=(512,512))
end
    mp4(anim,datadir("../films/soft_vision/real_seed_$(my_seed)/NRI_$(my_seed).mp4"))
    mp4(anim,datadir("../films/soft_vision/real_seed_$(my_seed)/NRI_$(my_seed)_fast.mp4"),fps=35)
# movie zoomed
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
locs = [(100,30)]
for (locx,locy) in locs
    z = @elapsed anim_zoom = @animate for tt in 1:length(saving_times)
        zoom_quiver(thetas_saves[:,:,tt],model,lattice,locx,locy,15,size=(1024,1024))
    end
    prinz(z)
    mp4(anim_zoom,datadir("../films/soft_vision/real_seed_$(my_seed)/NRI_$(my_seed)_Zoom($locx,$locy)_slow.mp4"),fps=15)
end
zoom_quiver(thetas_saves[:,:,100],model,lattice,locx,locy,15,size=(1024,1024))

## To have an idea of tmax
include(srcdir("../parameters.jl"));

tmax = 1E3
times = logspace(1,tmax,30)

R = 2
magnetisations = zeros(length(times),R)
corr_functions = zeros(Int(L/2),length(times),R)
corr_lengths   = zeros(length(times),R)
ns             = zeros(length(times),R)
z = @elapsed for r in 1:R
    println("Simulation $r/$R started.")
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
    token = 1
    while model.t < tmax # run the actual simulation
        update!(thetas,model,lattice)
        if model.t ≥ times[token] # once you hit a saving time, save quantites of interest
            println("t=$(round(model.t,digits=1)), $(round(model.t/tmax*100,digits=1))% achieved.")
            magnetisations[token,r]   = OP(thetas)[1] # polar Order Parameter
            corr_functions[:,token,r] = corr(thetas,model,lattice) # correlation function C(r)
            corr_lengths[token,r]     = corr_length(corr_functions[:,token,r]) # correlation length computed from C(r)
            ns[token,r] = number_defects(thetas,model,lattice)
            token = min(token+1,length(times)) # next saving_time
        end
    end
    p=plot_thetas(thetas,model,lattice,defects=false)
    display(p)
end
prinz(4/2*z)

# Average over the (possible) different realisations. Try with R=1 first, then increase if needed.
corr_functions_avg = mean(corr_functions,dims=3)
magnetisations_avg = mean(magnetisations,dims=2)
corrlength_avg = mean(corr_lengths,dims=2)
ns_avg = mean(ns,dims=2)
plot(times,ns_avg,axis=:log,m=true,xlabel="t",ylabel="n(t)")
    plot!(times,5E2 ./ (times).^0.66,c=:black)


## Test dt
using Random
Random.seed!(1)
include(srcdir("../parameters.jl"));
    model = SoftVisionXY(params)
    lattice = TriangularLattice(L)
    thetas = init_thetas(model,lattice,params_init=params_init)
z = @elapsed update!(thetas,model,lattice,tmax=100)
plot_thetas(thetas,model,lattice)


## Efficiency
using BenchmarkTools
model = SoftVisionXY(params)
model = XY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice)
@btime update!(thetas,model,lattice) # 14 ms versus 5 ms pour XY
