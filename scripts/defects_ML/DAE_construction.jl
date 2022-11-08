using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


using Flux, Augmentor, Random, Statistics, CUDA
using Flux:params, onehotbatch, crossentropy, logitcrossentropy, onecold, throttle, mse, trainmode!, testmode!
using Distributions, JLD2,BSON

#= --------------------- Idea ---------------------
One of the various uses of AutoEncoders is to denoise
samples. The NN basically learns to associate the noisy
sample to the "closest" to the "truth" it has learnt
(including from noisy samples only !)
Constraints : zero-mean noise, so that the leaning procedure
is unbaised.
=#

## Create base data set (without augmentation, just the different µ)
# include(srcdir("../parameters.jl"));
# dµ = pi/32
# mus = Float32.(round.(collect(0:dµ:2pi-dµ),digits=2))
# base_dataset = zeros(Float32,2WINDOW+1,2WINDOW+1,length(mus))
# model = XY(params)
# lattice = SquareLattice(W21,periodic=false)
# for i in each(mus)
#     params_init["type1defect"] = mus[i];
#     base_dataset[:,:,i] = init_thetas(model,lattice,params_init=params_init)
# end
# using JLD2
# jldsave("data/for_ML/base_dataset_µN12.jld2";base_dataset,mus,dµ,WINDOW)
# ind = rand(1:length(mus))
    # p=plot_thetas(base_dataset[:,:,ind],model,lattice)
    # display_quiver!(p,base_dataset[:,:,ind],WINDOW)

## Some functions
function random_flip(image,xpu;seuil=0.35)
    pi32 = Float32(π)
    return image .+ xpu(pi32*rand(Bernoulli(seuil), size(image)))
end

function proba_flip(epoch,epoch0=300,pmax=0.25;slope = 5E-4)
    linear_increase_at_from_e0 = max(0,slope*(epoch-epoch0))
    ceiled = min(linear_increase_at_from_e0,pmax)
    return ceiled
end
#
# function infer_mu(thetas::Array{T,4},q;window=WINDOW)::T where T<:AbstractFloat
#     thetas_reshaped = reshape(thetas,(2*window+1,2*window+1,1,1))
#     return infer_mu(thetas_reshaped,q,window=window)
# end
# function infer_mu(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
#     L = size(thetas,1)
#     @assert L == 2window+1
#     muss = zeros(size(thetas))
#     # tmp = Complex(0)
#     for j in 2:L-1, i in 2:L-1
#         muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
#         # i<->j irrelevant because i,j and j,i have the same weight for "mean" operation
#         # tmp += exp(im*muss[i,j] - 0.25sqrt((i-window)^2 + (j-window)^2))
#     end
#     # muss[window,window] = 0
#     moyenne = angle(mean(exp.(im*muss[2:L-1,2:L-1])))
#     if     abs(q) == 1   correction = pi - 0.33228605
#     elseif abs(q) == 1/2 correction = 0.1
#     end
#     #= tbh, I don't know where the shifts might come from,
#     In the beggining, I thought maybe from the spin at the center of the defect,
#     where theta should not be defined. But if one changes mean->sum and adds the condition
#     "if (i == j == window)", the result only becomes weirder... =#
#     # return mod.(muss,2pi)
#     return mod(moyenne .+ correction,2π)
# end
#
# # test infer_mu()
# inferred = zeros(length(mus))
#     for ind in 1:64
#         inferred[ind] = infer_mu(base_dataset[:,:,ind],1)
#     end
#     plot(mus,inferred,m=true)
#     plot!(x->x)

## Define Neural Network
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

#
function DenseResAutoEncoder(input_dim,hidden_dim,latent_dim)
    return Chain(
        Dense(input_dim,hidden_dim,relu),
        Dropout(0.5),
        SkipConnection(Dense(hidden_dim,hidden_dim,relu),+),
        SkipConnection(Dense(hidden_dim,hidden_dim,relu),+),
        Dense(hidden_dim,latent_dim,relu),
        Dropout(0.5),
        # . . . . . . . .

        Dense(latent_dim,hidden_dim,relu),
        Dropout(0.5),
        SkipConnection(Dense(hidden_dim,hidden_dim,relu),+),
        SkipConnection(Dense(hidden_dim,hidden_dim,relu),+),
        Dense(hidden_dim,input_dim,relu)
    )
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
export_to_gpu = true
    if export_to_gpu && CUDA.has_cuda()
        @info "Training on GPU"
        xpu = gpu
    else
        @info "Training on CPU"
        xpu = cpu
    end

## Data
W21 = 2WINDOW+1
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µN1.jld2")
include(srcdir("../parameters.jl"));

## Declare NN, Loss and Optimiser
L1norm(x) = sum(abs, x);
L1penalty() = sum(L1norm,Flux.params(NN))
L2norm(x) = sum(abs2, x);
L2penalty() = sum(L2norm,Flux.params(NN))
loss_pen(X, y) = mse(NN(X), y) + 5E-3*L2penalty()
loss(X, y) = mse(NN(X), y)
# progress = () -> @show(loss(Xtrain, Ytrain)) # callback to show loss
# evalcb = throttle(progress, 1)
NN = 0
    NN = ConvAutoEncoder(10) |> xpu
    opt = Adam(1E-3)
    # NN = DenseAutoEncoder(W21*W21,100,16) |> xpu

## Training
epochs = Int(6E3)
trainL = zeros(epochs)
trainLpen = zeros(epochs)

multi_fact = 5
params["symmetry"] = "polar"
model = LangevinXY(params)
lattice = TriangularLattice(W21,periodic=false)

X_noisy = similar(repeat(base_dataset,outer=[1,1,multi_fact]))
z = @elapsed for e in 1:epochs
    shuffled_dataset = repeat(base_dataset,outer=[1,1,multi_fact])[:,:,shuffle(1:end)]
    e0_noise = 1500 ; pmax = 0.3 ; slope = pmax/abs(epochs-e0_noise)*2
    seuil_flip = proba_flip(e,e0_noise,pmax,slope=slope)
    pi32 = Float32(π)
    for i in 1:size(shuffled_dataset,3)
        # Rotate
        if rand() < 0. degree = rand([0,10,20,30,350,340,330])
        else degree = rand(0:10:350)
        end
        ppl = Rotate(degree) |> Resize(W21,W21)
        tmp = augment(shuffled_dataset[:,:,i],ppl)
        tmp .+= Float32(deg2rad(degree))

        # Flip
        tmp .+= pi32*rand(Bernoulli(seuil_flip), size(tmp))

        # Thermal Noise (useless if relaxation at T>0)
        # tmp += 0.2*rand()*randn(size(tmp))

        # Relax
        # tmp = shuffled_dataset[:,:,i]
        trelax = .1 ; update!(tmp,model,lattice,trelax)

        X_noisy[:,:,i] = tmp
        # X_noisy[:,:,i] = shuffled_dataset[:,:,i] + Float32.(0.2*rand()*randn(W21,W21))
    end
    pi232 = Float32(2pi)
    X_noisy = mod.(X_noisy,pi232)
    # if needed, reshape to feed the conv/Dense layer
    X_noisy_reshaped = xpu(reshape(X_noisy,(W21,W21,1,:)))
    dataset_reshaped = xpu(reshape(shuffled_dataset,(W21,W21,1,:)))

    Flux.train!(loss_pen, Flux.params(NN),[(X_noisy_reshaped,dataset_reshaped)], opt)
    trainL[e] = loss(X_noisy_reshaped,dataset_reshaped)
    trainLpen[e] = loss_pen(X_noisy_reshaped,dataset_reshaped)

    if isinteger(e/10)
        println("e = $e/$epochs : train loss = $(round(trainL[e],digits=5))")
    end
end
prinz(z)

plot(legend=:bottomleft)
    plot!(1:epochs-1,trainL[1:end-1],axis=:log,lw=0.5)
    plot!(1:epochs-1,trainLpen[1:end-1],axis=:log,lw=0.5)
    plot!(1:epochs-1,(trainLpen - trainL)[1:end-1],axis=:log,lw=1)

# xx = reshape(base_dataset,(15,15,64,1))
# recon = NN(xx[:,:,1:1,1:1])[:,:]
# model = XY(params)
# lattice = TriangularLattice(15)
# plot_thetas(recon,model,lattice,defects=false)
#
# infer_mu(recon,1)
# NN = cpu(NN)

# comments = ["if trainL[e] < 0.5 opt.eta = 5E-4
#     elseif trainL[e] < 0.1 opt.eta = 2.5E-4
#     elseif trainL[e] < 0.02 opt.eta = 1E-4
#     end", "L1 1E-5 penalty, latent space dim = 10", "rotations in 0:10:350"]
using BSON
NN_saved = cpu(NN)
BSON.@save "DAE_negative12___07_11_2022.bson" NN_saved trainL trainLpen base_dataset epochs runtime=z
# @unpack NN, trainL, epochs = load(".jld2")

## END OF FILE
