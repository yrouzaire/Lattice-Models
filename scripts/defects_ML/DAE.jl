using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)

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

function random_flip(image,xpu;seuil=0.35)
    pi32 = Float32(π)
    return image .+ xpu(pi32*rand(Bernoulli(1-seuil), size(image)))
end

function add_noise(image,degree,seuil,xpu)
    image = image .+ Float32(deg2rad(degree)) # compensate for the rotation
    image = random_flip(image,xpu,seuil=seuil)
    return image
end

function get_proba_flip(epoch,epoch0=1000,pmax=0.3;slope = 5E-4)
    linear_increase_at_from_e0 = max(0,slope*(epoch-epoch0))
    ceiled = min(linear_increase_at_from_e0,pmax)
    return ceiled
end

function infer_mu(thetas::Matrix{T},q;window=WINDOW)::T where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    for j in 1:L, i in 1:L
        # muss[i,j] = thetas[i,j] - q*atan( (i-window) ,(j-window)) # i<->j doesn't change anything...
        muss[i,j] = thetas[i,j] - q*atan( (j-window) ,(i-window)) # i<->j doesn't change anything...
    end
    moyenne = angle(mean(exp.(im*muss)))
    # Correction des biais (aucune idée d'où ils sortent, ils ne sont pas vraiment constants non plus)
    # if q == +1   moyenne += 0.03 end
    # if q == -1/2 moyenne += 0.5 end
    # if q == +1/2 moyenne += 0.05 end
    return mod(moyenne,2π)
end

## Roadmap
# Generate sample defects (∀µ)
# Add "noise" (flip spins randomly) during training, more and more with epochs
# Test NN in denoising mode
# Use trivial routine to infer µ

## Define Neural Network
&
    function DenseAutoEncoder(input_dim,hidden_dim,latent_dim)
        encoder = Chain(
        Dense(input_dim,hidden_dim,relu),
        # Dropout(0.5),
        Dense(hidden_dim,latent_dim,relu)
        # Dropout(0.5)
        )
        decoder = Chain(
        Dense(latent_dim,hidden_dim,relu),
        # Dropout(0.5),
        Dense(hidden_dim,input_dim,relu),
        )
        autoencoder = Chain(encoder,decoder)
        return autoencoder
    end

# We define a reshape layer to use in our decoder
# Source :https://github.com/alecokas/flux-vae/blob/master/conv-vae/main.jl
    struct Reshape
        shape
    end
    Reshape(args...) = Reshape(args)
    (r::Reshape)(x) = reshape(x, r.shape)
    Flux.@functor Reshape ()

    function ConvAutoEncoder(latent_dim=16)
        output_conv_layer = (11,11,32)
        encoder = Chain(
            Conv((3, 3), 1=>16, relu),
            Conv((3, 3), 16=>32, relu),
                Flux.flatten,
            Dense(prod(output_conv_layer), 120, relu),
            Dense(120, latent_dim,relu)
            )
        decoder = Chain(
        Dense(latent_dim,120,relu),
        Dense(120,prod(output_conv_layer), relu),
        Reshape(output_conv_layer...,:),
        ConvTranspose((3,3),32=>16,relu),
        ConvTranspose((3,3),16=>1)
        )
        autoencoder = Chain(encoder,decoder)
        return autoencoder
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
opt = Adam()
loss(X, y) = mse(NN(X), y)
# progress = () -> @show(loss(Xtrain, Ytrain)) # callback to show loss
# evalcb = throttle(progress, 1)
NN = 0
    NN = ConvAutoEncoder(16) |> xpu
    # NN = DenseAutoEncoder(W21*W21,100,16) |> xpu

# NN(reshape(base_dataset,(15,15,1,63))[:,:,:,1:1])
## Training
base_dataset = xpu(base_dataset)
X_aug = similar(base_dataset)
drot = 10 ; rotation_angles = 0:drot:360-drot
epochs = Int(1E3)
trainL = zeros(epochs)
z = @elapsed for e in 1:epochs
    # Shuffle
    base_dataset[:,:,shuffle(1:end)]
    # Augmentation (Rotation)
    degree = 0#rand(rotation_angles)
    ppl_aug = Rotate(degree) |> Resize(W21,W21)
    augmentbatch!(X_aug,base_dataset,ppl_aug)
    # Add noise
    seuil = 0#get_proba_flip(e)
    X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
    # if needed, reshape to feed the conv/Dense layer
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    X_aug_reshaped = reshape(X_aug,(W21,W21,1,:))

    # X_noisy_reshaped = reshape(X_noisy,(W21*W21,:))
    # X_aug_reshaped = reshape(X_aug,(W21*W21,:))


    Flux.train!(loss, Flux.params(NN),[(X_noisy_reshaped,X_aug_reshaped)], opt)
    trainL[e] = loss(X_noisy_reshaped,X_aug_reshaped)

    println("e = $e/$epochs : train loss = $(round(trainL[e],digits=5))")
    # here, you may tune opt.η
    if trainL[e] < 0.5 opt.eta = 5E-4
    elseif trainL[e] < 0.1 opt.eta = 2.5E-4
    elseif trainL[e] < 0.02 opt.eta = 1E-4
    end
end
prinz(z)

plot!(1:10:epochs,trainL[1:10:end],axis=:log)


## Testing on noisy defects
model = XY(params)
lattice = TriangularLattice(W21)

degree = 0#rand(rotation_angles)
ppl_aug = Rotate(degree) |> Resize(W21,W21)
augmentbatch!(X_aug,base_dataset,ppl_aug)

ind = 60
    seuil = 0.0
    X_aug_reshaped = reshape(X_aug,(W21,W21,1,:))
    X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    X_noisy_reshaped += 0.2randn(size(X_noisy_reshaped))
    original = X_noisy_reshaped[:,:,1,ind]
        p1=plot_thetas(original,model,lattice,defects=false)
        display_quiver!(p1,original,WINDOW)
    recon = NN(X_noisy_reshaped[:,:,:,ind:ind])[:,:,1,1] .+ pi
    p2=plot_thetas(recon,model,lattice,defects=false)
    display_quiver!(p2,recon,WINDOW)
    plot(p1,p2,size=(485,800),layout=(2,1))

&
