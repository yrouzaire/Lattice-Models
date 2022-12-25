cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));

## Stability of µ over time
tmax = 10 ; every = 2 ; times = collect(0:every:tmax+every)
µ0s = [0,pi/2,pi,3pi/2]
µ0s = collect(0:pi/32:2pi)
reals = 10
visions = [0,0.05,0.1,0.2,0.4]
visions = [0.4]
mus = NaN*zeros(length(µ0s),length(visions),length(times)-1,reals)
z = @elapsed for i in each(µ0s)
    for j in each(visions)
        println("µ0 = ",round(µ0s[i],digits=1)," ; vision = ",visions[j])
        Threads.@threads for r in 1:reals
            params["L"] = 200 ; params["T"] = 0.1
            params["symmetry"] = "polar" ; params["rho"] = 1 ; params["vision"] = visions[j]
            params_init["init"] = "single" ; params_init["type1defect"] =  µ0s[i] # initial µ
            params_init["q"] = -1
            model = SoftVisionXY(params)
            lattice = TriangularLattice(L)
            thetas = init_thetas(model,lattice,params_init=params_init)

            try
                dft = DefectTracker(thetas,model,lattice,find_type=true)
                update_and_track!(thetas,model,lattice,dft,tmax,every,find_type=true)
                mus[i,j,:,r] = dft.defectsN[1].type[1:length(times)-1]
            catch ; end
        end
    end
end
prinz(z) # 5 minutes for a lot of curves, very fast
# @save "data/NRI/decay_scan_mu0.jld2" mus µ0s reals visions tmax every times runtime=z
# @load "data/NRI/decay_mu_someµ0_several_sigmas.jld2" mus µ0s reals visions tmax every times runtime
mus_avg = nanmean(mus,4)[:,:,:,1]
p=plot(ylims=(0,2pi),xlabel="t",ylabel="µ",legend=:outerright,size=(600,400))
    for j in 5#(each(visions))
    plot!([NaN,NaN],label="σ=$(visions[j])",c=j,rib=0)
    for i in 1:length(µ0s)
            plot!(times[1:end-1],mus_avg[i,j,:],c=j,lw=2)#,label="σ = $(visions[j])",rib=0)
            # plot!(visions[j]*times[1:end-1],mus_avg[i,j,:],c=j,lw=2)#,label="σ = $(visions[j])",rib=0)
        end
    end
    p
# savefig("plots/NRI/µt.png")
p=plot(ylims=(0,2pi),xlabel="σt",ylabel="µ")
    for i in [3,4]#:length(µ0s)
        for j in reverse(each(visions))
            plot!(visions[j]*times[1:35],mus_avg[i,j,1:35],c=j,lw=2)#,label="σ = $(visions[j])",rib=0)
        end
    end
    # plot!(0.1times[5:135],pi/2*exp.(-1.5(0.1times[5:135].-0.2)).+pi/2,c=:black,lw=2)
    # plot!(0.1times[5:135],pi/2*exp.(-1.5(0.1times[5:135].-0.2)).+pi/2,c=:black,lw=2)
    plot!(t->2atan(tanh(-(t-0.2)))+pi,c=:black,lw=2,line=:dash)
    p
# savefig("plots/NRI/µt_collapse.png")
p=plot(ylims=(0,2pi),xlabel="t",ylabel="µ")
    for i in 10:length(µ0s)-10
        for j in (each(visions))
            plot!(times[1:end-1],mus_avg[i,j,:],c=i,lw=2)#,label="σ = $(visions[j])",rib=0)
            # plot!(visions[j]*times[1:end-1],mus_avg[i,j,:],c=j,lw=2)#,label="σ = $(visions[j])",rib=0)
        end
    end
    plot!(times[1:end-1],mus_avg[49,1,:],c=:black,lw=2,line=:dash)#,label="σ = $(visions[j])",rib=0)
    plot!(times[1:end-1],mus_avg[17,1,:],c=:black,lw=2,line=:dashdot)#,label="σ = $(visions[j])",rib=0)
    p
# savefig("plots/NRI/µ_attractor.png")

#= The energy landscape seems to be something like -sin(µ)
Let's try and see what -σsin(µ) would give. =#
using DifferentialEquations
tspan = (0.0,10.0)
    f(u,p,t) = (0.4cos((u)))
    p=plot(legend=false,ylims=(0,2pi),ylabel="µ(t)")
    for u0 in 0:pi/16:2pi
        prob = ODEProblem(f,u0,tspan)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
        plot!(sol,lw=1)
        # cst = atanh(tan(µ0/2)/100)
        # plot!(t->2*atan(tanh(cst + t/2 )))
    end
    plot!(t->3pi/2,c=:black,lw=2,line=:dash)#,label="σ = $(visions[j])",rib=0)
    plot!(t->pi/2,c=:black,lw=2,line=:dashdot)#,label="σ = $(visions[j])",rib=0)
    p
# savefig("plots/NRI/traj_analytic_decayµ.png")
p=plot(legend=false,ylims=(0,2pi),ylabel="µ(t)",size=(400,400))
    tspan = (0.0,100.0)
    for sig in [0,0.05,0.1,0.2,0.4]
        f(u,p,t) = sig*cos(u)
        prob = ODEProblem(f,Float64(pi),tspan)
        sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
        plot!(sol,lw=1.5,col=sig,label="σ = $sig",rib=0)
    end
    # plot!(t->2atan(tanh(-(0.15*t)))+pi,c=:black,lw=2,line=:dash)
    p
# savefig("plots/NRI/traj_analytic_decayµ_sink.png")

## Movies
include(srcdir("../parameters.jl"));

    lattice = TriangularLattice(L)
    model = SoftVisionXY(params)
    thetas = init_thetas(model,lattice,params_init=params_init)
        plot_thetas(thetas,model,lattice)
    saving_times = 0:1:200 ; transients = saving_times[end]/2
    z = @elapsed anim = movies(thetas,model,lattice,defects=false,saving_times=saving_times,transients=transients)
    prinz(z)
    mp4(anim,datadir("../films/NRI/stability_defects/NRI_µ$(params["type1defect"]).mp4"))

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