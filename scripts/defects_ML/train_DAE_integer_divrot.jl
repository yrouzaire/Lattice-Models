using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

using Flux, Augmentor, Random, Statistics, CUDA
using Flux:params, onehotbatch, crossentropy, logitcrossentropy, onecold, throttle, mse, trainmode!, testmode!
using Distributions,BSON

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

## Define Neural Network
# We define a reshape layer to use in our decoder
# Source :https://github.com/alecokas/flux-vae/blob/master/conv-vae/main.jl
struct Reshape2
    shape
end
    Reshape2(args...) = Reshape2(args)
    (r::Reshape2)(x) = reshape(x, r.shape)
    Flux.@functor Reshape2 ()

    function ConvAE_divrot(latent_dim=16)
        output_conv_layer = (9,9,32)
        encoder = Chain(
            Conv((3, 3), 3=>16, relu),
            Conv((3, 3), 16=>32, relu),
            Conv((3, 3), 32=>32, relu),
                Flux.flatten,
            Dense(prod(output_conv_layer), 120, relu),
            Dropout(0.5),
            Dense(120,60,relu),
            Dropout(0.5),
            Dense(60,latent_dim,relu)
            )
        decoder = Chain(
        Dense(latent_dim,60,relu),
        Dense(60,120,relu),
        Dropout(0.5),
        Dense(120,prod(output_conv_layer), relu),
        Dropout(0.5),
        Reshape2(output_conv_layer...,:),
        ConvTranspose((3,3),32=>32,relu),
        ConvTranspose((3,3),32=>16,relu),
        ConvTranspose((3,3),16=>1)
        )
        autoencoder = Chain(encoder,decoder)
        return autoencoder
    end

    function ConvResAE_divrot(latent_dim=16)
        output_conv_layer = (9,9,32)
        return Chain(
        Conv((3, 3), 3=>16, relu),
        SkipConnection( # beggining of outer skip
            Chain(Conv((3, 3), 16=>32, relu),
            Conv((3, 3), 32=>32, relu),
            Flux.flatten,
            Dense(prod(output_conv_layer), 120, relu),
            SkipConnection(Chain( # beggining of inner skip
                Dense(120,60,relu),
                Dense(60,latent_dim,relu),
                # Latent Space #
                Dense(latent_dim,60,relu),
                Dense(60,120,relu)),+), # end of inner skip
            Dense(120,prod(output_conv_layer), relu),
            Reshape2(output_conv_layer...,:),
            ConvTranspose((3,3),32=>32,relu), # end of outer skip
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
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP1.jld2")
include(srcdir("../parameters.jl"));

## Declare NN, Loss and Optimiser
L1norm(x) = sum(abs, x); L1penalty() = sum(L1norm,Flux.params(NN))
L2norm(x) = sum(abs2, x); L2penalty() = sum(L2norm,Flux.params(NN))
loss_pen(X, y) = mse(NN(X), y) + 5E-3*L2penalty()
loss(X, y) = mse(NN(X), y)
# progress = () -> @show(loss(Xtrain, Ytrain)) # callback to show loss
# evalcb = throttle(progress, 1)
## Training
dim_latent_space = 3
NN = 0
    NN = ConvAE_divrot(dim_latent_space) |> xpu
    opt = Adam(1E-3)
    epochs = 2000

trainL = zeros(epochs)
    trainLpen = zeros(epochs)
    testL = zeros(epochs)

multi_fact = 10
X_noisy = similar(repeat(base_dataset,outer=[1,1,multi_fact]))
Ntrain = round(Int,0.8*size(X_noisy,3))
z = @elapsed for e in 1:epochs
    trainmode!(NN)
    shuffled_dataset = repeat(base_dataset,outer=[1,1,multi_fact])[:,:,shuffle(1:end)]

    for i in 1:size(shuffled_dataset,3) # Both Train and Validation Sets
        X_noisy[:,:,i] = shuffled_dataset[:,:,i] + Float32.(0.3 *rand()*randn(W21,W21))
    end
    pi232 = Float32(2pi)
    X_noisy = mod.(X_noisy,pi232)

    # Provide div and rot to help the NN
    X_aug_divrot = provide_div_rot(X_noisy)

    # if needed, reshape to feed the conv/Dense layer
    Xtrain = xpu(X_aug_divrot[:,:,:,1:Ntrain])
    Xtest = xpu(X_aug_divrot[:,:,:,1+Ntrain:end])
    Ytrain = xpu(reshape(shuffled_dataset[:,:,1:Ntrain],(W21,W21,1,:)))
    Ytest  = xpu(reshape(shuffled_dataset[:,:,1+Ntrain:end],(W21,W21,1,:)))

    Flux.train!(loss, Flux.params(NN),[(Xtrain,Ytrain)], opt)
    trainL[e] = loss(Xtrain,Ytrain)
    # trainLpen[e] = loss_pen(Xtrain,Ytrain)
    testmode!(NN)
    testL[e] =  loss(Xtest,Ytest)

    if isinteger(e/50)
        println("e = $e/$epochs : train loss = $(round(trainL[e],digits=5))")
    end
end
    prinz(z)

# plot(legend=:bottomleft)
    # plot!(1:epochs-1,trainLpen[1:end-1],axis=:log,lw=0.5,label="MSE + L2")
    # plot!(1:epochs-1,(trainLpen - trainL)[1:end-1],axis=:log,lw=1,label="L2")
plot!((1:epochs-1).+2000,testL[1:end-1],axis=:log,lw=0.5,label="Test")
    plot!((1:epochs-1).+2000,trainL[1:end-1],axis=:log,lw=0.5,label="MSE")

# comments = ["", "L1 1E-5 penalty, latent space dim = 10", "rotations in 0:10:350"]
using BSON
DAE = cpu(NN)
BSON.@save "NeuralNets/DAE_positive1___09_11_2022.bson" DAE trainL testL trainLpen base_dataset epochs runtime=z

## END OF FILE
