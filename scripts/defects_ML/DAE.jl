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
    return image .+ xpu(pi32*rand(Bernoulli(seuil), size(image)))
end

function add_noise(image,degree,seuil,xpu)
    image = image .+ Float32(deg2rad(degree)) # compensate for the rotation
    image = random_flip(image,xpu,seuil=seuil)
    return image
end

function proba_flip(epoch,epoch0=300,pmax=0.25;slope = 5E-4)
    linear_increase_at_from_e0 = max(0,slope*(epoch-epoch0))
    ceiled = min(linear_increase_at_from_e0,pmax)
    return ceiled
end

function rotation_angles(e,e0;slope=1E-3)
    if e < e0 return 0:0
    else
        drots = [180,90,45,40,35,30,25,20,10,5]
        drot = drots[min(length(drots),Int(floor(slope*(e-e0)))+1)]
        return 0:drot:360-drot
    end
end
function rotation_angles1(e,e0;slope=5E-3)
    drot = 10
    angle_max = min(drot*round(Int,max(0,slope*(e-e0)+1)),180)
    return (-angle_max:drot:angle_max) .+ 180
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
function DenseDropoutAutoEncoder(input_dim,hidden_dim,latent_dim)
    encoder = Chain(
    Dense(input_dim,hidden_dim,relu),
    Dropout(0.5),
    SkipConnection(Dense(hidden_dim,hidden_dim,relu),+),
    Dense(hidden_dim,latent_dim,relu),
    Dropout(0.5)
    )
    decoder = Chain(
    Dense(latent_dim,hidden_dim,relu),
    Dropout(0.5),
    SkipConnection(Dense(hidden_dim,hidden_dim,relu),+),
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
            Dense(120,60,relu),
            Dense(60,latent_dim,relu)
            )
        decoder = Chain(
        Dense(latent_dim,60,relu),
        Dense(60,120,relu),
        Dense(120,prod(output_conv_layer), relu),
        Reshape(output_conv_layer...,:),
        ConvTranspose((3,3),32=>16,relu),
        ConvTranspose((3,3),16=>1)
        )
        autoencoder = Chain(encoder,decoder)
        return autoencoder
    end

    function ConvResAutoEncoder0(latent_dim=16)
        output_conv_layer = (11,11,32)
        return Chain(

        Conv((3, 3), 1=>16, relu),
        SkipConnection( # beggining of outer skip
            Chain(Conv((3, 3), 16=>32, relu),
            Flux.flatten,
            Dense(prod(output_conv_layer), 120, relu),
            SkipConnection(Chain( # beggining of inner skip
                Dense(120,60,relu),
                Dense(60,latent_dim,relu),
                # Latent Space #
                Dense(latent_dim,60,relu),
                Dense(60,120,relu)),+), # end of inner skip
            Dense(120,prod(output_conv_layer), relu),
            Reshape(output_conv_layer...,:),
            ConvTranspose((3,3),32=>16,relu)),+), # end of outer skip
        ConvTranspose((3,3),16=>1)
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
sqnorm(x) = sum(abs2, x);
penalty() = sum(sqnorm,Flux.params(NN))
loss(X, y) = mse(NN(X), y) #+ 3E-5*penalty()
# progress = () -> @show(loss(Xtrain, Ytrain)) # callback to show loss
# evalcb = throttle(progress, 1)
NN = 0
    NN = ConvResAutoEncoder0(16) |> xpu
    opt = Adam()
    # NN = DenseAutoEncoder(W21*W21,100,16) |> xpu

# NN(reshape(base_dataset,(15,15,1,63))[:,:,:,1:1])

## Training
base_dataset = xpu(base_dataset)
epochs = Int(3E4)
trainL = zeros(epochs)

multi_fact = 10
X_aug = similar(repeat(base_dataset,outer=[1,1,multi_fact]))
z = @elapsed for e in 1:epochs
    # Shuffle
    shuffled_dataset = repeat(base_dataset[:,:,shuffle(1:end)],outer=[1,1,multi_fact])
    # Augmentation (Rotation)
    e0_angle = 400 ;
    degree = rand(rotation_angles1(e,e0_angle))
    ppl_aug = Rotate(degree) |> Resize(W21,W21)
    augmentbatch!(X_aug,shuffled_dataset,ppl_aug)

    # Add noise
    e0_noise = 3000 ; pmax = 0.3 ; slope = 1.5E-5
    seuil = proba_flip(e,e0_noise,pmax,slope=slope)
    X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
    # if needed, reshape to feed the conv/Dense layer
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    dataset_reshaped = reshape(shuffled_dataset,(W21,W21,1,:))

    # X_noisy_reshaped = reshape(X_noisy,(W21*W21,:))
    # dataset_reshaped = reshape(shuffled_dataset,(W21*W21,:))


    Flux.train!(loss, Flux.params(NN),[(X_noisy_reshaped,dataset_reshaped)], opt)
    trainL[e] = loss(X_noisy_reshaped,dataset_reshaped)

    if isinteger(e/50)
        println("e = $e/$epochs : train loss = $(round(trainL[e],digits=5))")
    end
    # here, you may tune opt.η
    if trainL[e] < 0.5 opt.eta = 5E-4
    elseif trainL[e] < 0.1 opt.eta = 2.5E-4
    elseif trainL[e] < 0.02 opt.eta = 1E-4
    end
end
prinz(z)
# using JLD2
# JLD2.@save "DAE_rotation+flip_positive12.jld2" NN trainL base_dataset

plot(legend=:bottomleft)
plot!(1:1:epochs,trainL[1:1:end],axis=:log,label="multi_fact = $multi_fact")

## Testing on noisy defects
model = XY(params)
lattice = TriangularLattice(W21)

ind = rand(1:63)
    drot = 1
    degree = rand(0:drot:360-drot)
    X_aug = similar(base_dataset)
    ppl_aug = Rotate(degree) |> Resize(W21,W21)
    augmentbatch!(X_aug,base_dataset,ppl_aug)
    seuil = 0.35
    X_aug_reshaped = reshape(X_aug,(W21,W21,1,:))
    # X_aug += 0.2randn(size(X_aug))
    X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    X_noisy_reshaped += 0.2randn(size(X_noisy_reshaped))
    thetas = base_dataset[:,:,ind]
        p0=plot_thetas(thetas,model,lattice,defects=false)
        display_quiver!(p0,thetas,WINDOW)
    original = X_noisy_reshaped[:,:,1,ind]
        p1=plot_thetas(original,model,lattice,defects=false)
        display_quiver!(p1,original,WINDOW)
    recon = NN(X_noisy_reshaped[:,:,:,ind:ind])[:,:,1,1]
    p2=plot_thetas(recon,model,lattice,defects=false)
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485*3,400),layout=(1,3))
infer_mu(mod.(recon,2pi),1/2)
mus[ind]

## Test on vivo defects
params["symmetry"] = "nematic" ; params["rho"] = 0.95 ; params["A"] = 2
model = SPP(params)
lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice,tmax=2000)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)

~,thetas_zoom= zoom(thetas,lattice,51,72)
zoom_quiver(thetas_zoom,model,lattice,8,8)
recon = reshape(NN(reshape(thetas_zoom,(15,15,1,1))),(15,15))
zoom_quiver(recon,model,lattice,8,8)
