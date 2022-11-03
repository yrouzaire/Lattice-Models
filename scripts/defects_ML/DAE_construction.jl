using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("../parameters.jl")); # do it before loading Flux (conflict with Flux.params)
include(srcdir("LatticeModels.jl"))

using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()


using Flux, Augmentor, Random, Statistics, CUDA
using Flux:params, onehotbatch, crossentropy, logitcrossentropy, onecold, throttle, mse, trainmode!, testmode!
using Distributions, JLD2

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

# function add_noise(image,degree,seuil,xpu)
#     image = image .+ Float32(deg2rad(degree)) # compensate for the rotation
#     image = random_flip(image,xpu,seuil=seuil)
#     return image
# end

function proba_flip(epoch,epoch0=300,pmax=0.25;slope = 5E-4)
    linear_increase_at_from_e0 = max(0,slope*(epoch-epoch0))
    ceiled = min(linear_increase_at_from_e0,pmax)
    return ceiled
end

# function rotation_angles(e,e0;slope=1E-3)
#     if e < e0 return 0:0
#     else
#         drots = [180,90,45,40,35,30,25,20,10,5]
#         drot = drots[min(length(drots),Int(floor(slope*(e-e0)))+1)]
#         return 0:drot:360-drot
#     end
# end
# function rotation_angles1(e,e0;slope=5E-3)
#     drot = 10
#     angle_max = min(drot*round(Int,max(0,slope*(e-e0)+1)),180)
#     return (-angle_max:drot:angle_max) .+ 180
# end

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
    # correction (might come from the non-correction of the rotation)
    # moyenne += 0.1
    return mod(moyenne,2π)
end

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
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP12.jld2")

## Declare NN, Loss and Optimiser
sqnorm(x) = sum(abs, x);
penalty() = sum(sqnorm,Flux.params(NN))
loss_pen(X, y) = mse(NN(X), y) + 1E-5*penalty()
loss(X, y) = mse(NN(X), y)
# progress = () -> @show(loss(Xtrain, Ytrain)) # callback to show loss
# evalcb = throttle(progress, 1)
NN = 0
    NN = ConvResAutoEncoder0(10) |> xpu
    opt = Adam(1E-3)
    # NN = DenseAutoEncoder(W21*W21,100,16) |> xpu

## Training
epochs = Int(2E4)
trainL = zeros(epochs)

multi_fact = 5
X_noisy = similar(repeat(base_dataset,outer=[1,1,multi_fact]))
z = @elapsed for e in 1:epochs
    shuffled_dataset = repeat(base_dataset,outer=[1,1,multi_fact])[:,:,shuffle(1:end)]
    e0_noise = 2000 ; pmax = 0.3 ; slope = pmax/abs(epochs-e0_noise)*2
    seuil_flip = 0#proba_flip(e,e0_noise,pmax,slope=slope)
    pi32 = Float32(π)
    for i in 1:size(shuffled_dataset,3)
        # Rotate
        degree = rand(0:10:350)
        ppl = Rotate(degree) |> Resize(W21,W21)
        tmp = augment(shuffled_dataset[:,:,i],ppl)
        tmp .+= Float32(deg2rad(degree))

        # Flip
        tmp .+= pi32*rand(Bernoulli(seuil_flip), size(tmp))
        X_noisy[:,:,i] = tmp
    end
    pi232 = Float32(2pi)
    X_noisy = mod.(X_noisy,pi232)
    # if needed, reshape to feed the conv/Dense layer
    X_noisy_reshaped = xpu(reshape(X_noisy,(W21,W21,1,:)))
    dataset_reshaped = xpu(reshape(shuffled_dataset,(W21,W21,1,:)))

    Flux.train!(loss_pen, Flux.params(NN),[(X_noisy_reshaped,dataset_reshaped)], opt)
    trainL[e] = loss(X_noisy_reshaped,dataset_reshaped)

    if isinteger(e/100)
        println("e = $e/$epochs : train loss = $(round(trainL[e],digits=5))")
    end
    # here, you may tune opt.η
    if trainL[e] < 0.5 opt.eta = 5E-4
    elseif trainL[e] < 0.1 opt.eta = 2.5E-4
    elseif trainL[e] < 0.02 opt.eta = 1E-4
    end
end
prinz(z)

comments = ["if trainL[e] < 0.5 opt.eta = 5E-4
    elseif trainL[e] < 0.1 opt.eta = 2.5E-4
    elseif trainL[e] < 0.02 opt.eta = 1E-4
    end","e0_noise = 300 ; pmax = 0.3 ; slope = pmax/(epochs-e0_noise)*1.5"]
JLD2.@save "DAE_positive12_02_11_2022.jld2" NN trainL base_dataset comments epochs runtime=z
@unpack NN, trainL, epochs = load("DAE_rotation+flip_positive12.jld2")

# plot(legend=:bottomleft)
plot!(1:epochs-1,trainL[1:end-1],axis=:log,lw=0.1)


## Testing (for CNN) on noisy defects
model = XY(params)
W21 = 2WINDOW+1
lattice = TriangularLattice(W21)

NN = cpu(NN)
ind = rand(1:63)
    degree = rand(0:10:350)
    ppl = Rotate(degree) |> Resize(W21,W21)
    tmp = augment(base_dataset[:,:,ind],ppl)
    tmp .+= Float32(deg2rad(degree))
    seuil_flip = 0
    tmp .+= Float32(pi)*rand(Bernoulli(seuil_flip), size(tmp))
    X_noisy = tmp
    # X_noisy += 0.2randn((W21,W21))
    pi232 = Float32(2pi)
    X_noisy = mod.(X_noisy,pi232)
    X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
    thetas = base_dataset[:,:,ind]
        p0=plot_thetas(thetas,model,lattice,defects=false)
        display_quiver!(p0,thetas,WINDOW)
    p1=plot_thetas(X_noisy,model,lattice,defects=false)
        display_quiver!(p1,X_noisy,WINDOW)
    recon = NN(X_noisy_reshaped)[:,:,1,1]
    p2=plot_thetas(recon,model,lattice,defects=false)
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485*3,400),layout=(1,3))
mus[ind]
infer_mu(mod.(recon,2pi),1/2)
&

## Estimating the error margin
R = 10
flip_strength = [0.1,0.2,0.3]
errors = zeros(length(mus),length(flip_strength),R)

drot = 10 ; ppl_aug = Rotate(0:drot:360-drot) |> Resize(W21,W21)
z = @elapsed for i in each(mus)
    original = base_dataset[:,:,i]
    for j in each(flip_strength) , r in 1:R
        X_aug = similar(original)
        X_aug = augment(original,ppl_aug)
        seuil = flip_strength[j]
        X_aug_reshaped = reshape(X_aug,(W21,W21,1,:))
        X_aug += 0.2randn(size(X_aug))
        X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
        X_noisy_reshaped = reshape(X_noisy,(W21,W21,1,:))
        recon = NN(X_noisy_reshaped[:,:,:,1:1])[:,:,1,1]
        errors[i,j,r] = abs(infer_mu(mod.(recon,2pi),1/2) .+ 0.1 - mus[i])
    end
end
prinz(z)
errors_avg = mean(errors,dims=3)[:,:,1]
p=plot(legend=:topleft)
    for j in each(flip_strength)
        plot!(mus[4:end],errors_avg[4:end,j],label="Proba Flip = $(flip_strength[j])")
    end
    p

## Testing for Dense AutoEncoder
ind = rand(1:63)
    drot = 1
    degree = rand(0:drot:360-drot)
    X_aug = similar(base_dataset)
    ppl_aug = Rotate(degree) |> Resize(W21,W21)
    augmentbatch!(X_aug,base_dataset,ppl_aug)
    seuil = 0.
    # X_aug += 0.2randn(size(X_aug))
    X_noisy = mapslices(image->add_noise(image,degree,seuil,xpu),X_aug,dims=(1,2))
    # X_noisy += 0.2randn(size(X_noisy))
    thetas = base_dataset[:,:,ind]
        p0=plot_thetas(thetas,model,lattice,defects=false)
        display_quiver!(p0,thetas,WINDOW)
    original = X_noisy[:,:,ind]
        p1=plot_thetas(original,model,lattice,defects=false)
        display_quiver!(p1,original,WINDOW)
    recon = reshape(NN(vec(X_noisy[:,:,ind])),(W21,W21))
    p2=plot_thetas(recon,model,lattice,defects=false)
    display_quiver!(p2,recon,WINDOW)
    plot(p0,p1,p2,size=(485*3,400),layout=(1,3))
mus[ind]
infer_mu(mod.(recon,2pi),1/2)

## Test on vivo defects :
params["symmetry"] = "nematic" ; params["rho"] = 1 ; params["A"] = 0
lattice = TriangularLattice(L)
model = LangevinXY(params)

thetas = init_thetas(model,lattice,params_init=params_init)
update!(thetas,model,lattice) # calentamiento
update!(thetas,model,lattice,tmax=20)  # updates until time = t
    p = plot_thetas(thetas,model,lattice,defects=false)

    dft = DefectTracker(thetas,model,lattice,find_type=true)
    histogram(last_types(dft),bins=50,normalize=true)
update_and_track!(thetas,model,lattice,dft,2000,100,find_type=true)
number_active_defects(dft)
filter(isnan,last_types(dft))

~,thetas_zoom= zoom(thetas,lattice,114,120)
zoom_quiver(thetas_zoom,model,lattice,8,8)
recon = reshape(DAE_positive12(reshape(thetas_zoom,(15,15,1,1))),(15,15))
zoom_quiver(recon,model,lattice,8,8)
infer_mu(recon,1/2)
