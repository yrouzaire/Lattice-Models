cd("D:/Documents/Research/projects/LatticeModels")
    using DrWatson ; @quickactivate "LatticeModels"
    include(srcdir("LatticeModels.jl"))
    using Plots,ColorSchemes,LaTeXStrings
    pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
    include(srcdir("../parameters.jl"));
&
CHARGE = 1
    @unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP1.jld2")
    base_datasetP = base_dataset
    @unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µN1.jld2")
    base_datasetN = base_dataset;

## Visually check the dataset
ind = 17#rand(1:length(mus))
    p=plot_thetas(base_datasetP[:,:,ind],model,lattice)
    display_quiver!(p,base_datasetP[:,:,ind],WINDOW)
    title!("µ = $(mus[ind])")

ind = 33#rand(1:length(mus))
    p=plot_thetas(base_datasetN[:,:,ind],model,lattice)
    display_quiver!(p,base_datasetN[:,:,ind],WINDOW)
    title!("µ = $(mus[ind])")

## Test noiseless
inferredP  = zeros(length(mus))
    inferredN = zeros(length(mus))
    decay = true
    for ind in 1:64
        inferredP[ind] = infer_mu(base_datasetP[:,:,ind],q=CHARGE,decay=decay)
        inferredN[ind] = infer_mu(base_datasetN[:,:,ind],q=-CHARGE,decay=decay)
    end
    p1 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noiseless")
    plot!(mus,inferredP,m=true,c=:red,label="q=1")
    plot!(mus,inferredN,m=true,c=:green,label="q=-1")
    plot!(x->x,c=:black)

## Test with noise
inferredP  = zeros(length(mus))
    inferredN = zeros(length(mus))
    noise = 0.3randn(W21,W21,64)
    decay = false
    for ind in 1:64
        inferredP[ind] = infer_mu(base_datasetP[:,:,ind] + noise[:,:,ind],q=CHARGE,decay=decay)
        inferredN[ind] = infer_mu(base_datasetN[:,:,ind]+ noise[:,:,ind],q=-CHARGE,decay=decay)
    end
    p2 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noisy (without decay)")
    plot!(mus,inferredP,line=false,m=true,c=:red,label="q=1")
    plot!(mus,inferredN,line=false,m=true,c=:green,label="q=-1")
    plot!(x->x,c=:black)

decay = true
    for ind in 1:64
        inferredP[ind] = infer_mu(base_datasetP[:,:,ind] + noise[:,:,ind],q=CHARGE,decay=decay)
        inferredN[ind] = infer_mu(base_datasetN[:,:,ind] + noise[:,:,ind],q=-CHARGE,decay=decay)
    end
    p3 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noisy (with decay for q=1)")
    plot!(mus,inferredP,line=false,m=true,c=:red,label="q=1")
    plot!(mus,inferredN,line=false,m=true,c=:green,label="q=-1")
    plot!(x->x,c=:black)

plot(p1,p2,p3,size=(1200,400),layout=(1,3))
# savefig("plots/procedure_infer_mu_polar.png")
