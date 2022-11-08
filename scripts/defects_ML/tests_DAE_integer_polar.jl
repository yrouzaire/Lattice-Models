using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));
params["symmetry"] = "polar"
CHARGE = 1
model = XY(params)
lattice = TriangularLattice(W21,periodic=false)


using CUDA, Flux, BSON
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP1.jld2")
@unpack NN, trainL, epochs = load("DAE_positive1___07_11_2022_fonctionnel.bson")
BSON.@load "DAE_positive1___07_11_2022_fonctionnel.bson" NN_saved
NN_test = cpu(NN_saved)
# comments

using Augmentor

## First Test : reconstruct randomly generated noisy, flipped and rotated in vitro defect.
params["symmetry"] = "polar"
model = XY(params)
lattice = TriangularLattice(W21,periodic=false)
ind = rand(1:63)
    tmp = base_dataset[:,:,ind]
    trelax = .5 ; update!(tmp,model,lattice,trelax)
    X_noisy = tmp
    pi232 = Float32(2pi)
        X_noisy = mod.(X_noisy,pi232)
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    thetas = base_dataset[:,:,ind]
        p0=plot_thetas(thetas,model,lattice,defects=false,title="µ = $(mus[ind])")
        display_quiver!(p0,thetas,WINDOW)
    p1=plot_thetas(X_noisy,model,lattice,defects=false)
        display_quiver!(p1,X_noisy,WINDOW)
    recon = NN_test(X_noisy_reshaped)[:,:,1,1]
    p2=plot_thetas(recon,model,lattice,defects=false,title="µ = $(round(infer_mu(recon,q=CHARGE),digits=2))")
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485,400*3),layout=(3,1))

mus[ind]
round(infer_mu(mod.(recon,2pi),q=CHARGE),digits=2)
round(infer_mu(recon,q=CHARGE),digits=2)
heatmap((thetas-recon)')
thetas_modified = thetas .+ 0.2randn(size(thetas))*rand(Bernoulli(0.1), size(thetas))
    pp = plot_thetas(thetas_modified,model,lattice,defects=false)
    display_quiver!(pp,thetas_modified,WINDOW)

round(infer_mu(mod.(thetas,2pi),q=CHARGE),digits=2)
round(infer_mu(mod.(thetas_modified,2pi),q=CHARGE),digits=2)

## Second Test : check the reconstruction is correct for all mus
R = 100
mus_infered = zeros(length(mus),R)

z = @elapsed for i in each(mus)
    original = base_dataset[:,:,i]
    for r in 1:R
        tmp = original
        trelax = .1 ; update!(tmp,model,lattice,trelax)
        tmp = mod.(tmp,2pi)

        X_noisy_reshaped = reshape(tmp,(W21,W21,1,:))
        recon = NN_test(X_noisy_reshaped)[:,:,1,1]
        mus_infered[i,r] = infer_mu_decay(tmp,q=CHARGE)
    end
end
plot(xlabel="True µ",ylabel="Inferred µ")
    scatter!(mus,mus_infered[1:end,:],c=:black,ms=1)
    plot!(x->x-0.3,c=:red,lw=0.5)
    plot!(x->x+0.3,c=:red,lw=0.5)

## Third Test : Errors when inferring true µ = 0
R = 1000
flip_strength = [0.1,0.2,0.3]
mus_infered = zeros(length(flip_strength),R)
original = base_dataset[:,:,1] # µ = 0


z = @elapsed for j in each(flip_strength) , r in 1:R
    tmp = original
    trelax = .1 ; update!(tmp,model,lattice,trelax)
    tmp = mod.(tmp,2pi)
    X_noisy_reshaped = reshape(tmp,(W21,W21,1,:))
    recon = NN_test(X_noisy_reshaped)[:,:,1,1]
    mus_infered[j,r] = infer_mu(recon,q=CHARGE)
end
histogram(mod.(mus_infered[2,:].-pi,2pi).+pi,bins=50,normalize=true)
mean(mod.(mus_infered[1,:].+pi,2pi).-pi)
std(mod.(mus_infered[1,:].+pi,2pi).-pi)

## Forth Test : check the constance over time with in vitro defects
p=plot(ylims=(0,2pi))
    for initmu in 1:6 ,  r in 1:1
        params["L"] = 200 ; params["T"] = 0.2
        params["symmetry"] = "polar" ; params["rho"] = 1 ; params["A"] = 0
        params_init["init"] = "single" ; params_init["type1defect"] =  initmu # initial µ
        params_init["q"] = CHARGE
        model = MonteCarloXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        # p = plot_thetas(thetas,model,lattice,defects=false)

        dft = DefectTracker(thetas,model,lattice,find_type=true)
        update_and_track!(thetas,model,lattice,dft,250,10,find_type=true)
        plot!(dft.defectsP[1].type,m=false,line=true,c=initmu)
    end
    p
&

## Fifth Test : In vivo defects, LangevinXY, measure P(µ,t)
include(srcdir("../parameters.jl"));
params["symmetry"] = "polar" ; params["rho"] = 1 ; params["A"] = 0
params["L"] = 200 ; params["T"] = 0.2
lattice = TriangularLattice(L)
model = MonteCarloXY(params)

thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,tmax=100)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)
thetasagg = zeros(Float32,2 .*size(thetas))
    thetasagg[1:L,1:L] = thetasagg[1:L,1+L:end] = thetasagg[1+L:end,1:L] = thetasagg[1+L:end,1+L:end] = thetas
    latticeagg = TriangularLattice(2L)
update!(thetasagg,model,latticeagg,tmax=500)  # updates until time = t
    p = plot_thetas(thetasagg,model,latticeagg,defects=false)

dft = DefectTracker(thetasagg,model,latticeagg,find_type=true)
    plot(xlims=(0,2pi))
    histogram!(last_types(dft),bins=28,normalize=true)

loc = dft.defectsP[1].pos[1]
    zoom_quiver(thetasagg,model,latticeagg,loc...)
~,thetas_zoom= zoom(thetasagg,latticeagg,loc...)
    recon = reshape(NN_test(reshape(thetas_zoom,(15,15,1,1))),(15,15))
    zoom_quiver(recon,model,TriangularLattice(W21),8,8)
infer_mu(recon,q=CHARGE)