cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

#= This file aims at testing the correctness of MSD(dft) on known cases, namely XY/ForcedXY pairs and isolated defects.
XY Single : seems OK, starts with ~t increase and then sort of plateaus, weird
Forced XY Single : seems OK, starts with ~t increase and then sort of plateaus, weird

XY Pair : seems OK, starts with ~t increase and then sort of plateaus, weird
=#

include(srcdir("../parameters.jl"));
    lattice = SquareLattice(L)
    R = 10
    dfts = Vector{Union{Missing,DefectTracker}}(missing,R)
    tmax, every = 250,5

z = @elapsed Threads.@threads for r in 1:R
    println()
    println("r = $r/$R")
    model = XY(params)
    thetas = init_thetas(model,lattice,params_init=params_init)
    dft = DefectTracker(thetas,model,lattice)
    update_and_track!(thetas,model,lattice,dft,tmax,every)
    dfts[r] = dft
end
prinz(z)
msd = MSD(dfts,lattice)[1]
    plot(msd[3:end],axis=:log)
plot!(x->10/100*x)
ind += 1
    R = interdefect_distance(dfts[ind],dfts[ind].defectsP[1],dfts[ind].defectsN[1],lattice)
    plot!(R)
