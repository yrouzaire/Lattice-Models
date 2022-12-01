cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings,SpecialFunctions,LambertW
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

filename = datadir("first_simus_NRI.jld2")
@unpack runtimes,scanned_params, params, comments, R, times_log, times_lin  = load(filename)
@unpack Ts,visions = scanned_params
comments
L = params["L"]
lattice = TriangularLattice(L)

## Coarsening dynamics
@unpack polar_orders, nematic_orders, Cs, xis, ns = load(filename)
polar_orders_avg = nanmean(polar_orders,4)[:,:,:,1]
nematic_orders_avg = nanmean(nematic_orders,4)[:,:,:,1]
ns_avg  = nanmean(ns,4)[:,:,:,1]
xis_avg = nanmean(xis,4)[:,:,:,1]
Cs_avg  = nanmean(Cs,5)[:,:,:,:,1]
histogram(runtimes/3600/24,bins=40)


ls = [:solid,:dash]
# Order parameters (t)
p=plot(xlabel="t",ylabel="OP",legend=:topleft,axis=:log,size=(450,400))
    for i in each(Ts)
        for j in 1:1#each(visions)
            plot!(times_log,polar_orders_avg[i,j,:],c=i,line=:solid,label="T = $(Ts[i])",rib=0)
            plot!(times_log,nematic_orders_avg[i,j,:],c=i,line=:dash)
        end
    end
    plot!(times_log[2:end],5E-2sqrt.(times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Correlation length (t)
p=plot(xlabel="t",ylabel="ξ/L",legend=false,axis=:log,size=(450,400))
    for i in each(Ts)
        for j in 1:1#each(visions)
            plot!(times_log,xis_avg[i,j,:]/L,c=i,line=ls[j],label="T = $(Ts[i])",rib=0)
        end
    end
    plot!(times_log[2:end],2E-2sqrt.(times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Number of defects (t)
p=plot(xlabel="t",ylabel="n",legend=:bottomleft,axis=:log,size=(450,400))
    for i in each(Ts)
        for j in 1:2#each(visions)
            plot!(times_log,ns_avg[i,j,:].+1,c=i,line=ls[j],label="T = $(Ts[i])",rib=0)
        end
    end
    plot!(times_log[2:end],1E4 ./ (times_log[2:end]./log.(5times_log[2:end])),c=:black)
    # plot!(times_log[2:end],2E3 ./ (times_log[2:end]./log.(5times_log[2:end]).^2),c=:black)
    p

# Correlation function over time
temp = 3
    p=plot(xlabel="r",ylabel="C(r,t) for T=$(Ts[temp])",legend=:bottomleft,axis=:log,size=(450,400))
    for i in 27:1:31
        plot!(1:length(Cs_avg[temp,:,i]),remove_negative(Cs_avg[temp,:,i]),line=:solid,label="t = $(times_log[i])",rib=0)
    end
    p

# Correlation function at final time for different T
p=plot(xlabel="r",ylabel="C(r,∞)",legend=:bottomleft,axis=:log,size=(450,400))
    for i in 1:length(Ts)
        rr = 1:length(Cs_avg[i,:,end])
        plot!(rr,remove_negative(Cs_avg[i,:,end]),c=i,line=:solid,label="T = $(Ts[i])",rib=0)
        plot!(rr[10:end],.97rr[10:end] .^(-Ts[i]/pi/1.5),c=i,line=:dash)
    end
    p

## DefectTracker
@unpack dftss = load(filename)

p=plot(legend=:topleft)
    for i in each(Ts)
        msd = MSD(dftss[i,:],lattice)[1]
        plot!(msd[3:end],axis=:log,label="T = $(Ts[i])",rib=0)
    end
    p
    plot!(x->15(x/log(10x)),c=:black)

p=plot(legend=:topleft)
    for i in each(Ts)
        R = mean_distance_to_annihilator(dftss[i,:],lattice)
        plot!(R[2:end],axis=:log,label="T = $(Ts[i])",rib=0)
    end
    p
    plot!(x->2exp(lambertw(4pi*x)/2),c=:black)
    # plot!(x->3sqrt(x),c=:black)
    # plot!(x->10sqrt(x/log(10x)),c=:black)

## Thetas saved
@unpack thetas_saves = load(filename)
model = XY(params)
thetass = thetas_saves[1,1,24,:,:,1]
    plot_thetas(thetass,model,lattice,defects=false)
zoom_quiver(thetass,model,lattice,115,110,10)

anim = @animate for i in 1:size(thetas_saves,3)
    plot_thetas(thetas_saves[1,1,i,:,:,1],model,lattice,defects=false,size=(512,512))
end
mp4(anim,"films\\NRI_test00.mp4",fps=5)