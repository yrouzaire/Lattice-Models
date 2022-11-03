using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


using Flux, Augmentor, Random, Statistics, CUDA
using Flux:params, onehotbatch, crossentropy, logitcrossentropy, onecold, throttle, mse, trainmode!, testmode!
using Distributions

#= --------------------- Idea ---------------------
One of the various uses of AutoEncoders is to denoise
samples. The NN basically learns to associate the noisy
sample to the "closest" to the "truth" it has learnt
(including from noisy samples only !)
Constraints : zero-mean noise, so that the leaning procedure
is unbaised.
=#

## Roadmap
# Method that returns many channels with rotated images
# Add "noise" (flip spins randomly) during training, more and more with epochs
# Test NN in denoising mode
# Use trivial routine to infer µ

function image2batch_rotation(image;deg=5)
    degrees = 0:deg:360-deg
    batch = zeros(Float32,(size(image)...,length(degrees)))
    for i in each(degrees)
        degree = degrees[i]
        ppl = Rotate(degree) |> Resize(size(image,1),size(image,2))
        batch[:,:,i] = augment(image,ppl) .+ Float32(deg2rad(degree))
    end
    return batch[:,:,shuffle(1:end)]
end

function create_batch(images;deg=5,permutation_mus)
    a,b,c = size(images)
    degrees = 0:deg:360-deg
    batch = zeros(Float32,(a,b,length(degrees),c)) # space x, space y, rotations, µ
    for i in 1:c
        batch[:,:,:,i] = image2batch_rotation(images[:,:,i])
    end
    return batch[:,:,:,permutation_mus] # shuffles the µ
end

function random_flip(image,xpu;seuil=0.35)
    pi32 = Float32(π)
    return image .+ xpu(pi32*rand(Bernoulli(seuil), size(image)))
end

function proba_flip(epoch,epoch0=300,pmax=0.25;slope = 5E-4)
    linear_increase_at_from_e0 = max(0,slope*(epoch-epoch0))
    ceiled = min(linear_increase_at_from_e0,pmax)
    return ceiled
end

function infer_mu(thetas::Array{T,4},q;window=WINDOW)::T where T<:AbstractFloat
    thetas_reshaped = reshape(thetas,(2*window+1,2*window+1,1,1))
    return infer_mu(thetas_reshaped,q,window=window)
end
function infer_mu(thetas::Matrix{T},q;window=WINDOW)::T where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    for j in 1:L, i in 1:L
        muss[i,j] = thetas[i,j] - q*atan( (i-window) ,(j-window))
        # i<->j irrelevant because i,j and j,i has the same weight for "mean" operation
    end
    moyenne = angle(mean(exp.(im*muss)))
    return mod(moyenne,2π)
end



## Define Neural Network
&
# We define a reshape layer to use in our decoder
# Source :https://github.com/alecokas/flux-vae/blob/master/conv-vae/main.jl
    struct Reshape
        shape
    end
    Reshape(args...) = Reshape(args)
    (r::Reshape)(x) = reshape(x, r.shape)
    Flux.@functor Reshape ()
    function CNN(nb_classes=64)
        output_conv_layer = (9,9,9) # prod = 729
        return Chain(
        Conv((3, 3), 72=>36, relu),
        Conv((3, 3), 36=>18, relu),
        Conv((3, 3), 18=>9, relu),
        Flux.flatten,
        Dense(prod(output_conv_layer), 120, relu),
        Dense(120,60,relu),
        Dense(60,nb_classes),
        softmax
        )
    end
## Device CPU/GPU
export_to_gpu = false
    if export_to_gpu && CUDA.has_cuda()
        @info "Training on GPU"
        xpu = gpu
    else
        @info "Training on CPU"
        xpu = cpu
    end

## Data
W21 = 2WINDOW+1
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP12.jld2")

## Declare NN, Loss and Optimiser
# sqnorm(x) = sum(abs2, x);
# penalty() = sum(sqnorm,Flux.params(NN))
loss(X, y) = crossentropy(NN(X), y) #+ 3E-5*penalty()
# progress = () -> @show(loss(Xtrain, Ytrain)) # callback to show loss
# evalcb = throttle(progress, 1)
NN = 0
    NN = CNN() |> xpu
    opt = Adam(1E-3)
    # NN = DenseAutoEncoder(W21*W21,100,16) |> xpu

# NN(reshape(base_dataset,(15,15,1,63))[:,:,:,1:1])

## Training
base_dataset = xpu(base_dataset)
epochs = Int(1E2)
trainL = zeros(epochs)

drot = 5 ; rotations = 0:drot:360-drot
X_fliped = zeros(Float32,(W21,W21,length(rotations),length(mus)))
z = @elapsed for e in 1:epochs
    # Shuffle
    perm1 = randperm(length(mus))
    X = create_batch(base_dataset,deg=drot,permutation_mus=perm1)
    Y = onehotbatch(mus[perm1],mus)

    # Augmentation (Flips)
    e0_noise = 1500 ; pmax = 0.25 ; slope = 1E-4
    seuil = proba_flip(e,e0_noise,pmax,slope=slope)
    X_fliped = mapslices(img->random_flip(img,xpu,seuil=seuil),X,dims=(1,2))

    Flux.train!(loss, Flux.params(NN),[(X_fliped,Y)], opt)
    trainL[e] = loss(X_fliped,Y)

    if isinteger(e)
        println("e = $e/$epochs : train loss = $(round(trainL[e],digits=5))")
    end
    # here, you may tune opt.η
    if trainL[e] < 5 opt.eta = 5E-7
    elseif trainL[e] < 0.1 opt.eta = 2.5E-4
    elseif trainL[e] < 0.02 opt.eta = 1E-4
    end
end
prinz(z)
using JLD2
JLD2.@save "DAE_divrot_positive12.jld2" NN trainL base_dataset

plot(legend=false)#:bottomleft)
plot!(1:100:epochs,trainL[1:100:end],axis=:log,label="multi_fact = $multi_fact")
plot(smooth(trainL[1:1:end];over = 100))

## Testing (for CNN) on noisy defects
model = XY(params)
W21 = 2WINDOW+1
lattice = TriangularLattice(W21)

ind = rand(1:63)
    drot = 15
    degree = rand(0:drot:360-drot)
    ppl_aug = Rotate(degree) |> Resize(W21,W21)
    X_aug = similar(base_dataset)
    augmentbatch!(X_aug,base_dataset,ppl_aug)
    seuil = 0.
    X_aug_reshaped = reshape(X_aug,(W21,W21,1,:))
    # X_aug += 0.2randn(size(X_aug))
    X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    # X_noisy_reshaped += 0.2randn(size(X_noisy_reshaped))

    X_noisy_divrot = zeros(Float32,W21,W21,3,1)
    X_noisy_divrot[:,:,1,1] = X_noisy[:,:,1] # ch1 = noisy thetas
    divmat, rotmat = get_div_rot(X_noisy[:,:,1],lattice)
    X_noisy_divrot[:,:,2,1] = divmat # ch2 = div
    X_noisy_divrot[:,:,3,1] = rotmat # ch3 = rot

    thetas = base_dataset[:,:,ind]
        p0=plot_thetas(thetas,model,lattice,defects=false)
        display_quiver!(p0,thetas,WINDOW)
    original = X_noisy_reshaped[:,:,1,ind]
        p1=plot_thetas(original,model,lattice,defects=false)
        display_quiver!(p1,original,WINDOW)
    recon = NN(X_noisy_divrot[:,:,:,1:1])[:,:,1,1]
    p2=plot_thetas(recon,model,lattice,defects=false)
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485*3,400),layout=(1,3))
mus[ind]
infer_mu(mod.(recon,2pi),1/2)

## Test on vivo defects : XY
params["symmetry"] = "nematic" ; params["rho"] = 1 ; params["A"] = 0
lattice = TriangularLattice(L)
model = MonteCarloXY(params)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,1000)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)
model.t

dft = DefectTracker(thetas,model,lattice,find_type=true)
update_and_track!(thetas,model,lattice,dft,2000,100,find_type=true)
histogram(last_types(dft),bins=50,normalize=true)
dft

~,thetas_zoom= zoom(thetas,lattice,114,120)
zoom_quiver(thetas_zoom,model,lattice,8,8)
recon = reshape(DAE_positive12(reshape(thetas_zoom,(15,15,1,1))),(15,15))
zoom_quiver(recon,model,lattice,8,8)
infer_mu(recon,1/2)
