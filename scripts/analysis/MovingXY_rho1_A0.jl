cd("D:/Documents/Research/projects/LatticeModels")
using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


#= Important Comments : =#
filename = datadir("polarMovXY_rho1_A0.jld2")
## Coarsening dynamics
@unpack runtimes, polar_orders, nematic_orders, Cs, xis, ns, params, comments, R, times_log, times_log, Ts  = load(filename)
polar_orders_avg = nanmean(polar_orders,3)[:,:,1]
nematic_orders_avg = nanmean(nematic_orders,3)[:,:,1]
ns_avg = nanmean(ns,3)[:,:,1]
xis_avg = nanmean(xis,3)[:,:,1]
Cs_avg = nanmean(Cs,4)[:,:,:,1]
L = params["L"]
histogram(runtimes/3600,bins=R)

# Order parameters (t)
p=plot(xlabel="t",ylabel="OP",legend=:topleft,axis=:log,size=(450,400))
    for i in each(Ts)
        plot!(times_log,polar_orders_avg[i,:],c=i,line=:solid,label="T = $(Ts[i])",rib=0)
        plot!(times_log,nematic_orders_avg[i,:],c=i,line=:dash)
    end
    plot!(times_log[2:end],5E-2sqrt.(times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Correlation length (t)
p=plot(xlabel="t",ylabel="ξ/L",legend=:topleft,axis=:log,size=(450,400))
    for i in each(Ts)
        plot!(times_log,xis_avg[i,:]/L,c=i,line=:solid,label="T = $(Ts[i])",rib=0)
    end
    plot!(times_log[2:end],2E-2sqrt.(times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Number of defects (t)
p=plot(xlabel="t",ylabel="n/L²",legend=:bottomleft,axis=:log,size=(450,400))
    for i in each(Ts)
        plot!(times_log,ns_avg[i,:].+1,c=i,line=:solid,label="T = $(Ts[i])",rib=0)
    end
    plot!(times_log[2:end],1E-2 ./ (times_log[2:end]./log.(10times_log[2:end])),c=:black)
    p

# Correlation function over time
temp = 3
    p=plot(xlabel="r",ylabel="C(r,t) for T=$(Ts[temp])",legend=:bottomleft,axis=:log,size=(450,400))
    for i in 20:2:30
        plot!(1:length(Cs_avg[temp,:,i]),remove_negative(Cs_avg[temp,:,i]),line=:solid,label="t = $(times_log[i])",rib=0)
    end
    p

# Correlation function at final time for different T
p=plot(xlabel="r",ylabel="C(r,∞)",legend=:bottomleft,axis=:log,size=(450,400))
    for i in 1:length(Ts)
        rr = 1:length(Cs_avg[i,:,end])
        plot!(rr,remove_negative(Cs_avg[i,:,end]),c=i,line=:solid,label="T = $(Ts[i])",rib=0)
        plot!(rr[10:end],.95rr[10:end] .^(-Ts[i]/pi/1.47),c=i,line=:dash)
    end
    p

## Thetas
@unpack thetas_saves = load(filename)
model = XY(params)
lattice = TriangularLattice(L)
plot_thetas((thetas_saves[1,25,:,:,1]),model,lattice,defects=false)
zoom_quiver((thetas_saves[1,25,:,:,1]),model,lattice,125,45,10)

## DefectTracker
@unpack dftss, times_lin = load(filename)
lattice = TriangularLattice(L)
dftss[1]

MSD_all, MSD_P, MSD_N = MSD(dftss[1,:],lattice)
plot(filter(!isnan,MSD_P[2:end]),axis=:log)
plot!(filter(!isnan,MSD_N[2:end]),axis=:log)
plot(filter(!isnan,MSD_all[2:end]),axis=:log)

indices = [] # indices of dft defined (is simulation not finished, dfts[i] == missing)
for i in 1:length(dftss[1,:])
    if !ismissing(dftss[1,:][i]) push!(indices,i) end
end
maxlength = maximum([maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)]) for dft in dftss[1,indices]])

dft = dftss[1,2]
    maximum([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)])

argmax([length(d.pos) for d in vcat(dft.defectsP,dft.defectsN)])

plot(remove_negative(square_displacement(dftss[1,2].defectsP[28],lattice)),axis=:log)
scatter(dftss[1,2].defectsP[28].pos)
scatter!(dftss[1,2].defectsP[28].pos[1],c=:red)
scatter!(dftss[1,2].defectsP[28].pos[end],c=:green)
dftss[1,2].defectsP[28].pos

plot(every:every:every*length(MSD_P)-1,MSD_P[2:end],axis=:log)
plot(every:every:every*length(MSD_N)-1,MSD_N[2:end],axis=:log)
plot(every:every:every*length(MSD_all)-1,MSD_all[2:end],axis=:log)
