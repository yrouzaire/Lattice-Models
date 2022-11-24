using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));
params["symmetry"] = "nematic"
CHARGE = 1/2
model = XY(params)
lattice = TriangularLattice(W21,periodic=false)

@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP12.jld2")
BSON.@load "NeuralNets/DAE_positive12___10_11_2022.bson" DAE
NN_test = cpu(DAE)
NN_test = cpu(NN)
# comments

using Augmentor

## First Test : reconstruct randomly generated noisy, flipped and rotated in vitro defect.
ind = rand(1:64)
    degree = rand(0:10:350)
        ppl = Rotate(degree) |> Resize(W21,W21)
        tmp = augment(base_dataset[:,:,ind],ppl)
        tmp .+= Float32(deg2rad(degree))
    seuil_flip = 0.
        tmp .+= Float32(pi)*rand(Bernoulli(seuil_flip), size(tmp))
    # trelax = .1 ; update!(tmp,model,lattice,trelax)
    X_noisy = tmp
    pi232 = Float32(2pi)
        X_noisy = mod.(X_noisy,pi232)
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    thetas = base_dataset[:,:,ind]
    p0=plot_thetas(thetas,model,lattice,defects=false,title="True µ = $(mus[ind])")
        display_quiver!(p0,thetas,WINDOW)
    p1=plot_thetas(X_noisy,model,lattice,defects=false,title="Inferred µ = $(round(infer_mu(X_noisy,q=CHARGE),digits=2))")
        display_quiver!(p1,X_noisy,WINDOW)
    recon = NN_test(provide_div_rot_muss(X_noisy_reshaped))[:,:,1,1]
    p2=plot_thetas(recon,model,lattice,defects=false,title="Inferred µ = $(round(infer_mu(recon,q=CHARGE),digits=2))")
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485,400*3),layout=(3,1))
    # plot(p0,p1,p2,size=(485*3,400),layout=(1,3))
    # savefig("plots\\DAE\\nematic\\recon.png")

## Second Test : check the reconstruction is correct for all mus
R = 100
    flip_strength = [0.]
    mus_infered_withDAE = zeros(length(mus),length(flip_strength),R)
    mus_infered_withoutDAE = zeros(length(mus),length(flip_strength),R)

    z = @elapsed for i in each(mus)
    original = base_dataset[:,:,i]
    for j in each(flip_strength) , r in 1:R
        degree = rand(0:10:350)
            ppl = Rotate(degree) |> Resize(W21,W21)
            tmp = augment(original,ppl)
            tmp .+= Float32(deg2rad(degree))

        tmp .+= Float32(pi)*rand(Bernoulli(flip_strength[j]), size(tmp))

        trelax = .1 ; update!(tmp,model,lattice,trelax)
        tmp = mod.(tmp,2pi)

        X_noisy_reshaped = reshape(tmp,(W21,W21,1,:))
        recon = NN_test(provide_div_rot_muss(X_noisy_reshaped))[:,:,1,1]
        mus_infered_withDAE[i,j,r] = infer_mu_decay(recon,q=CHARGE)
        mus_infered_withoutDAE[i,j,r] = infer_mu_decay(tmp,q=CHARGE)
    end
end
p1=plot(xlabel="True µ",ylabel="Inferred µ",title="Without DAE preprocess. ")
    scatter!(mus,mus_infered_withoutDAE[1:end,1,:],c=:black,ms=1)
    plot!(x->x-0.3,c=:red,lw=0.5)
    plot!(x->x+0.3,c=:red,lw=0.5)

    p2=plot(xlabel="True µ",ylabel="Inferred µ",title="With DAE preprocess. ")
    scatter!(mus,mus_infered_withDAE[1:end,1,:],c=:black,ms=1)
    plot!(x->x-0.3,c=:red,lw=0.5)
    plot!(x->x+0.3,c=:red,lw=0.5)

    plot(p1,p2,size=(800,400))
# savefig("plots/DAE/nematic/comparison_with_without_DAE_nematic.png")

## Third Test : Errors when inferring true µ = 0
R = 1000
flip_strength = [0.3]
mus_infered = zeros(length(flip_strength),R)
original = base_dataset[:,:,1] # µ = 0


z = @elapsed for j in each(flip_strength) , r in 1:R
    degree = rand(0:10:350)
        ppl = Rotate(degree) |> Resize(W21,W21)
        tmp = augment(original,ppl)
        tmp .+= Float32(deg2rad(degree))

    tmp .+= Float32(pi)*rand(Bernoulli(flip_strength[j]), size(tmp))

    trelax = .1 ; update!(tmp,model,lattice,trelax)
    tmp = mod.(tmp,2pi)
    X_noisy_reshaped = reshape(tmp,(W21,W21,1,:))
    recon = NN_test(X_noisy_reshaped)[:,:,1,1]
    mus_infered[j,r] = infer_mu(recon,q=CHARGE)
end
plot(xlabel="µ",ylabel=L"P(µ|µ_{true} = 0)")
    histogram!(mod.(mus_infered[1,:].+pi,2pi).-pi,bins=50,normalize=true)
    avg = round(mean(mod.(mus_infered[1,:].+pi,2pi).-pi),digits=2)
    stdd = round(std(mod.(mus_infered[1,:].+pi,2pi).-pi),digits=2)
    annotate!(-1.7,0.9,text("Mean = $avg",8))
    annotate!(-1.7,0.8,text("Std dev = $stdd",8))
# savefig("plots/DAE/nematic/error_at_µ0.png")

## Forth Test : check the constance over time with in vitro defects
p=plot(ylims=(0,2pi),xlabel=("t"),ylabel=("µ(t)"))
    for initmu in 1:5 ,  r in 1:1
        println("µ = $initmu , r = $r")
        params["L"] = 200 ; params["T"] = 0.2
        params["symmetry"] = "nematic" ; params["rho"] = 1 ; params["A"] = 0
        params_init["init"] = "single" ; params_init["type1defect"] =  initmu # initial µ
        params_init["q"] = CHARGE
        model = MonteCarloXY(params)
        lattice = TriangularLattice(L)
        thetas = init_thetas(model,lattice,params_init=params_init)
        thetas .+= Float32(pi)*rand(Bernoulli(0.2), size(thetas))
        # p = plot_thetas(thetas,model,lattice,defects=false)

        dft = DefectTracker(thetas,model,lattice,find_type=true)
        update_and_track!(thetas,model,lattice,dft,250,10,find_type=true)
        plot!(dft.defectsP[1].type,m=false,line=true,c=initmu)
    end
    p
&
# savefig("plots/DAE/nematic/constance_over_time_nematic.png")


## Fifth Test : In vivo defects, LangevinXY, measure P(µ,t)
include(srcdir("../parameters.jl"));
params["symmetry"] = "nematic" ; params["rho"] = 1 ; params["A"] = 0
params["L"] = 200 ; params["T"] = 0.2
lattice = TriangularLattice(L)
latticeagg = TriangularLattice(2L)

R = 10
dfts = Vector{DefectTracker}(undef,R)
z = @elapsed for r in 1:R
    println("r=$r/$R")
    model = MonteCarloXY(params)
    thetas = init_thetas(model,lattice,params_init=params_init)
    update!(thetas,model,lattice) # calentamiento
    update!(thetas,model,lattice,tmax=100)  # updates until time = t
        p = plot_thetas(thetas,model,lattice,defects=false)
    thetasagg = zeros(Float32,2 .*size(thetas))
        thetasagg[1:L,1:L] = thetasagg[1:L,1+L:end] = thetasagg[1+L:end,1:L] = thetasagg[1+L:end,1+L:end] = thetas
    update!(thetasagg,model,latticeagg,tmax=500)  # updates until time = t
        # p = plot_thetas(thetasagg,model,latticeagg,defects=false)

    dfts[r] = DefectTracker(thetasagg,model,latticeagg,find_type=true)
end
prinz(z)
last_types_vec = []
    for r in 1:R push!(last_types_vec,last_types(dfts[r])...) end
    plot(xlims=(0,2pi),xlabel="µ",ylabel="P(µ)")
    histogram!(last_types_vec,bins=50,normalize=true)
    hline!([1/2pi],c=:black)
# savefig("plots/DAE/nematic/Pmu_nematic.png")


loc = dft.defectsP[1].pos[1]
loc = (109,143)
    zoom_quiver(thetasagg,model,latticeagg,loc...)
~,thetas_zoom= zoom(thetasagg,latticeagg,loc...)
    recon = reshape(NN_test(reshape(thetas_zoom,(15,15,1,1))),(15,15))
    zoom_quiver(recon,model,TriangularLattice(W21),8,8)
infer_mu(thetas_zoom,q=CHARGE)
infer_mu(recon,q=CHARGE)


#
dft = DefectTracker(thetas,model,lattice,find_type=true)
ind = 8
    loc = dft.defectsP[ind].pos[1]
    loc =
    zoom_quiver(thetas,model,lattice,loc...)
    title!("Inferred µ = $(round(dft.defectsP[ind].type[1],digits=2))")
    # savefig("plots\\DAE\\nematic\\illustration_difficulty_inference.png")
