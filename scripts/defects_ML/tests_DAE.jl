using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));
model = XY(params)
lattice = TriangularLattice(W21)

@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP1.jld2")
@unpack NN, trainL, epochs = load("DAE_positive1___03_11_2022.jld2")
NN_test = cpu(NN)

## First Test : reconstruct randomly generated noisy, flipped and rotated in vitro defect.
ind = rand(1:63)
    degree = 0#rand(0:10:350)
    ppl = Rotate(degree) |> Resize(W21,W21)
    tmp = augment(base_dataset[:,:,ind],ppl)
    tmp .+= Float32(deg2rad(degree))
    seuil_flip = 0.0
    tmp .+= Float32(pi)*rand(Bernoulli(seuil_flip), size(tmp))
    X_noisy = tmp
    X_noisy += 0.2randn((W21,W21))
    pi232 = Float32(2pi)
    X_noisy = mod.(X_noisy,pi232)
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    thetas = base_dataset[:,:,ind]
        p0=plot_thetas(thetas,model,lattice,defects=false,title="µ = $(mus[ind])")
        display_quiver!(p0,thetas,WINDOW)
    p1=plot_thetas(X_noisy,model,lattice,defects=false)
        display_quiver!(p1,X_noisy,WINDOW)
    recon = NN_test(X_noisy_reshaped)[:,:,1,1]
    p2=plot_thetas(recon,model,lattice,defects=false,title="µ = $(round(infer_mu(mod.(recon,2pi),q=1),digits=2))")
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485,400*3),layout=(3,1))

mus[ind]
round(infer_mu(mod.(recon,2pi),q=1),digits=2)
round(infer_mu(recon,q=1),digits=2)
heatmap((thetas-recon)')
thetas_modified = thetas .+ 0.1randn(size(thetas))*rand(Bernoulli(0.1), size(thetas))
    pp = plot_thetas(thetas_modified,model,lattice,defects=false)
    display_quiver!(pp,thetas_modified,WINDOW)

round(infer_mu(mod.(thetas,2pi),q=1),digits=2)
round(infer_mu(mod.(thetas_modified,2pi),q=1),digits=2)
histogram(vec(infer_mu(mod.(thetas,2pi),q=1)),bins=50)
histogram!(vec(infer_mu(mod.(thetas_modified,2pi),q=1)),bins=50)
## Second Test : check the reconstruction is correct for all mus
plot_thetas(thetas,model,lattice,defects=false)


R = 100
flip_strength = [0.1,0.2,0.3]
flip_strength = [0.]
mus_infered = zeros(length(mus),length(flip_strength),R)

z = @elapsed for i in each(mus)
    original = base_dataset[:,:,i]
    for j in each(flip_strength) , r in 1:R
        # Rotate
        # degree = rand(0:10:350)
        # ppl = Rotate(degree) |> Resize(W21,W21)
        # tmp = augment(original,ppl)
        # tmp .+= Float32(deg2rad(degree))

        tmp = original
        # Flip
        # tmp .+= Float32(pi)*rand(Bernoulli(flip_strength[j]), size(tmp))

        # Thermal Noise
        tmp .+= Float32.(0.2*rand()*randn(size(tmp)))

        tmp = mod.(tmp,2pi)
        X_noisy_reshaped = reshape(tmp,(W21,W21,1,:))
        recon = NN_test(X_noisy_reshaped)[:,:,1,1]
        mus_infered[i,j,r] = infer_mu(mod.(recon,2pi),q=1)
    end
end
prinz(z)
plot(xlabel="True µ",ylabel="Inferred µ")
    scatter!(mus,mus_infered[1:end,1,:],c=:black,ms=1)
    plot!(x->x,c=:red)
    plot!(x->x-0.3,c=:red,lw=0.5)
    plot!(x->x+0.3,c=:red,lw=0.5)

## Third Test : Errors when inferring true µ = 0
R = 1000
flip_strength = [0.1,0.2,0.3]
mus_infered = zeros(length(flip_strength),R)
original = base_dataset[:,:,20] # µ = 0

z = @elapsed for j in each(flip_strength) , r in 1:R
        # Rotate
        degree = rand(0:10:350)
        ppl = Rotate(degree) |> Resize(W21,W21)
        tmp = augment(original,ppl)
        tmp .+= Float32(deg2rad(degree))

        # Flip
        tmp .+= Float32(pi)*rand(Bernoulli(flip_strength[j]), size(tmp))
        tmp = mod.(tmp,2pi)
        X_noisy_reshaped = reshape(tmp,(W21,W21,1,:))
        recon = NN_test(X_noisy_reshaped)[:,:,1,1]
        mus_infered[j,r] = infer_mu(mod.(recon,2pi),q=1/2)
end
histogram(mod.(mus_infered[3,:],2pi),bins=50,normalize=true)

## Forth Test : check the constance over time with in vitro defects
p=plot(ylims=(0,2pi))
    for initmu in 1:6 ,  r in 1:2
        params["symmetry"] = "polar" ; params["rho"] = 1 ; params["A"] = 0
        params_init["init"] = "single" ; params_init["type1defect"] =  initmu # initial µ
        # params_init["init"] = "pair" ; params_init["type2defect"] = [0,pi/2] # initial µ
        model = MonteCarloXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        thetas .+= Float32(pi)*rand(Bernoulli(0.), size(thetas))
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
lattice = TriangularLattice(L)
model = MonteCarloXY(params)

thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,tmax=200)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)
thetasagg = zeros(Float32,2 .*size(thetas))
    thetasagg[1:L,1:L] = thetasagg[1:L,1+L:end] = thetasagg[1+L:end,1:L] = thetasagg[1+L:end,1+L:end] = thetas
    latticeagg = TriangularLattice(2L)
update!(thetasagg,model,latticeagg,tmax=500)  # updates until time = t
    p = plot_thetas(thetasagg,model,latticeagg,defects=false)

dft = DefectTracker(thetasagg,model,latticeagg,find_type=true)
number_active_defects(dft)
plot(xlims=(0,2pi))
histogram!(last_types(dft),bins=50,normalize=true)

loc = dft.defectsP[1].pos[1]
~,thetas_zoom= zoom(thetas,lattice,loc...)
zoom_quiver(thetasagg,model,latticeagg,loc...)
recon = reshape(DAE_positive12(reshape(thetas_zoom,(15,15,1,1))),(15,15))
zoom_quiver(recon,model,lattice,8,8)
infer_mu(mod.(recon,2pi),q=1/2)
