cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
gr(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
include(srcdir("../parameters.jl"));


## Efficiency
using BenchmarkTools
model = SoftVisionXY(params)
model = XY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice)
@btime update!(thetas,model,lattice) # 14 ms versus 5 ms pour XY

## Movies
include(srcdir("../parameters.jl"));

model = SoftVisionXY(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
    plot_thetas(thetas,model,lattice)
saving_times = 0:1:200 ; transients = saving_times[end]/2
z = @elapsed anim = movies(thetas,model,lattice,defects=false,saving_times=saving_times,transients=transients)
prinz(z)
mp4(anim,datadir("../films/soft_vision/test3_NRI.mp4"))


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


## Analysis data cluster
@unpack scanned_params, params, runtimes, = load()
