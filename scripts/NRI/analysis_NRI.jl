cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings,SpecialFunctions,LambertW
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
ls = [:solid,:dash,:dot,:dashdot]

filename = datadir("first_simus_NRI.jld2")
@unpack runtimes,scanned_params, params, comments, R, times_log, times_lin  = load(filename)
@unpack Ts,visions = scanned_params
comments
L = params["L"]
lattice = TriangularLattice(L)

## Loading
@unpack polar_orders, nematic_orders, Cs, xis, ns = load(filename)
polar_orders_avg = nanmean(polar_orders,4)[:,:,:,1]
nematic_orders_avg = nanmean(nematic_orders,4)[:,:,:,1]
ns_avg  = nanmean(ns,4)[:,:,:,1]
xis_avg = nanmean(xis,4)[:,:,:,1]
Cs_avg  = nanmean(Cs,5)[:,:,:,:,1]
histogram([runtimes]/3600/24,bins=40)

## Phase Space (Final time and Movie over time)
# Number of defects
plot(xlabel="Vision",ylabel="T",colorbartitle=L"log_{10}(n+1)",size=(500,400))
    heatmap!(visions,Ts,log10.(ns_avg[:,:,end] .+ 1),c=cgrad([:blue,:green,:gold,:orangered,:red]),clims=(0,maximum(log10.(ns_avg .+ 1))))
anim = @animate for tt in 1:size(ns_avg,3)
    plot(xlabel="Vision",ylabel="T",colorbartitle=L"log_{10}(n+1)",size=(500,400))
    heatmap!(visions,Ts,log10.(ns_avg[:,:,tt] .+ 1),c=cgrad([:blue,:green,:gold,:orangered,:red]),clims=(0,maximum(log10.(ns_avg .+ 1))))
end
mp4(anim,"plots/NRI/phase_spaceTVision_over_time.mp4")

# Polar Order
plot(xlabel="Vision",ylabel="T",colorbartitle="P",size=(500,400))
    heatmap!(visions,Ts,polar_orders_avg[:,:,end],c=cgrad(reverse([:blue,:green,:gold,:orange,:red])),clims=(0,1))
anim = @animate for tt in 1:size(ns_avg,3)
    plot(xlabel="Vision",ylabel="T",colorbartitle="P",size=(500,400))
    heatmap!(visions,Ts,polar_orders_avg[:,:,tt],c=cgrad(reverse([:blue,:green,:gold,:orange,:red])),clims=(0,1))
end
mp4(anim,"plots/NRI/P_phase_spaceTVision_over_time.mp4")

## Coarsening dynamics
# Order parameters (t)
p=plot(xlabel="t",ylabel="OP",legend=false,axis=:log,size=(450,400))
    for i in 5#each(Ts)
        for j in each(visions)
            plot!(times_log,polar_orders_avg[i,j,:],c=i,line=:solid,label="Vision = $(visions[j])",rib=0)
            # plot!(times_log,nematic_orders_avg[i,j,:],c=i,line=:dash)
        end
    end
    plot!(times_log[2:end],5E-2sqrt.(times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Correlation length (t)
p=plot(xlabel="t",ylabel="ξ/L",legend=false,axis=:log,size=(450,400))
    for i in 2#each(Ts)
        for j in each(visions)
            plot!(times_log,xis_avg[i,j,:]/L,c=i,line=:solid,label="Vision = $(visions[j])",rib=0)
        end
    end
    plot!(times_log[2:end],2E-2sqrt.(times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Number of defects (t)
p=plot(xlabel="t",ylabel="n",legend=false,axis=:log,size=(450,400))
    for i in 7#each(Ts)
        for j in each(visions)
            plot!(times_log,ns_avg[i,j,:].+1,c=i,line=:solid,label="Vision = $(visions[j])",rib=0)
        end
    end
    plot!(times_log[2:end],1E3 ./ (times_log[2:end]./log.(5times_log[2:end])),c=:black)
    # plot!(times_log[2:end],2E3 ./ (times_log[2:end]./log.(5times_log[2:end]).^2),c=:black)
    p

# Correlation function over time
temp = 2 ; vis = 1
    p=plot(xlabel="r",ylabel="C(r,t) for T=$(Ts[temp])",legend=:outerright,axis=:log,size=(650,400))
    for i in 20:28
        plot!(1:length(Cs_avg[temp,vis,:,i]),remove_negative(Cs_avg[temp,vis,:,i]),line=:solid,label="t = $(times_log[i])",rib=0)
    end
    p

# Correlation function at final time for different T
p=plot(xlabel="r",ylabel="C(r,∞)",legend=:bottomleft,axis=:log,size=(450,400))
    for i in 1#:length(Ts)
        for j in each(visions)
            rr = 1:length(Cs_avg[i,j,:,end])
            plot!(rr,remove_negative(Cs_avg[i,j,:,end]),c=i,line=ls[j],label="Vision = $(visions[j])",rib=0)
            plot!(rr[10:end],.97rr[10:end] .^(-Ts[i]/pi/1.5),c=i,line=:dash)
        end
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
thetass = thetas_saves[2,5,25,:,:,1]
    plot_thetas(thetass,model,lattice,defects=false)
zoom_quiver(thetass,model,lattice,115,110,10)

anim = @animate for i in 1:size(thetas_saves,3)
    plot_thetas(thetas_saves[2,5,i,:,:,1],model,lattice,defects=false,size=(512,512))
end
mp4(anim,"films\\NRI_test00.mp4",fps=10)

## Separations µµ
filename = datadir("separations_µµ.jld2")
@load filename mus separationss sigmas r0 runtimes
rs = separationss
rs_avg = nanmean(rs,5)[:,:,:,:,1]/r0

plot(xlabel=L"µ_{+}",ylabel=L"µ_{-}",size=(470,400))
    heatmap!(mus,mus,rs_avg[:,:,2,end],c=cgrad([:blue,:deepskyblue2,:white,:red,:red2]),aspect_ratio=1,colorbartitle="R(t)/R0",clims=(minimum(0),1))
    # plot!(mus,5pi/2 .+ pi/2 .- mus)
savefig("plots/NRI/attraction_µµ_sigma0.3.png")
anim = @animate for tt in size(rs_avg,4)
    plot(xlabel=L"µ_{+}",ylabel=L"µ_{-}",size=(470,400))
        heatmap!(mus,mus,rs_avg[:,:,2,end],c=cgrad([:blue,:deepskyblue2,:white,:red,:red2]),aspect_ratio=1,colorbartitle="R(t)/R0",clims=(minimum(0),1))
end
mp4("plots/NRI/separations_µµ_sigma0.3.mp4")
