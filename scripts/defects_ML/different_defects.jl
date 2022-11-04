using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));

#= New Idea
In fact, there is a continuum of defects type.
θ = q*arctan(Δy/Δx) + µ, with 0 ≤ µ ≤ 2π.

Special cases of q>0 defects:
µ = 0 source
µ = π/2 counterclockwise
µ = π sink
µ = 3π/2 clockwise

Special cases of q<0 defects:
µ = 0 source
µ = π/2 counterclockwise
µ = π sink
µ = 3π/2 clockwise

This resolves the problem of having a defect
that looks like a counterclockwise AND a source
at the same time. It just means that 0 < µ < π/2.

=#

function reverse_spins(x::Matrix{T},p::Number)::Matrix{T} where T<:AbstractFloat
    N = length(x)
    indices = sample(1:N,round(Int,N*p),replace=false)
    pi32 = T(π)
    for ind in indices
        x[ind] += pi32
    end
    return mod.(x,2pi32)
end

## Now one has to retrieve µ from the WINDOW x WINDOW portion of the theta field
include(srcdir("../parameters.jl"));
# mu = 2pi*rand();
mu = 5;
    model = XY(params)
    lattice = SquareLattice(L)
    params_init["type1defect"] = mu;
    thetas = init_thetas(model,lattice,params_init=params_init)
    update!(thetas,model,lattice,5)
    window = WINDOW
    thetas_rotated = randomly_rotate(thetas)
    muss = infer_mu(zoom(thetas_rotated,lattice,spot_defects(thetas_rotated,model,
    lattice)[2][1][1:2]...,window)[2],q,window=window)
mu

function infer_mu(thetas::Matrix{T},q;window=WINDOW)::T where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    for j in 1:L, i in 1:L
        muss[i,j] = thetas[i,j] - q*atan( (i-window) ,(j-window)) # i<->j doesn't change anything...
        muss[i,j] = thetas[i,j] - q*atan( (j-window) ,(i-window)) # i<->j doesn't change anything...
    end
    moyenne = angle(mean(exp.(im*muss)))
    # Correction des biais (aucune idée d'où ils sortent, ils ne sont pas vraiment constants non plus)
    # if q == +1   moyenne += 0.03 end
    # if q == -1/2 moyenne += 0.5 end
    # if q == +1/2 moyenne += 0.05 end
    return mod(moyenne,2π)
end

## Get hands on on data augmentation (rotations + gaussian noise)
using Augmentor, Rotations, ImageTransformations

include(srcdir("../parameters.jl"));
mu = 0;
    model = XY(params)
    lattice = SquareLattice(L)
    params_init["type1defect"] = mu;
    thetas = init_thetas(model,lattice,params_init=params_init)

rotangle = 45
    pl = Rotate(rotangle) |> Resize(2WINDOW+1, 2WINDOW+1)
    thetas_rotated =  augment(thetas,pl) .+ Float32(deg2rad(rotangle)) + 0.25*randn(2WINDOW+1,2WINDOW+1)
    p=plot_thetas(thetas_rotated,model,lattice,defects=false)
    display_quiver!(p,thetas_rotated,WINDOW)

## Create base data set (without augmentation, just the different µ)
include(srcdir("../parameters.jl"));
dµ = pi/32
mus = Float32.(round.(collect(0:dµ:2pi-dµ),digits=2))
base_dataset = zeros(Float32,2WINDOW+1,2WINDOW+1,length(mus))
model = XY(params)
lattice = SquareLattice(L)
for i in each(mus)
    params_init["type1defect"] = mus[i];
    base_dataset[:,:,i] = init_thetas(model,lattice,params_init=params_init)
end
# using JLD2
# jldsave("data/for_ML/base_dataset_µP12.jld2";base_dataset,mus,dµ,WINDOW)
ind = rand(1:length(mus))
    p=plot_thetas(base_dataset[:,:,ind],model,lattice)
    display_quiver!(p,base_dataset[:,:,ind],WINDOW)

## Augmentation of the base_dataset for Dense NN
using Augmentor
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP12.jld2")
rot_angles = 0:15:285
    nb_noise = 200
    N = length(rot_angles)*nb_noise*length(mus)
    X = zeros(Float32,L*L,N)
    Y = zeros(Float32,N)
    token = 1
    for i in each(mus)
        mu = mus[i]
        for rotation_angle in rot_angles
            for n in 1:nb_noise
                ppl = Rotate(rotation_angle) |> Resize(L,L)
                # augmented_thetas = augment(base_dataset[:,:,i],ppl) .+ Float32(deg2rad(rotation_angle))
                augmented_thetas = augment(base_dataset[:,:,i],ppl) .+ Float32(deg2rad(rotation_angle)) + Float32(0.2)*rand(Float32)*randn(Float32,L,L)
                X[:,token] = vec(augmented_thetas)
                Y[token] = mu
                token += 1
            end
        end
    end
X
Y


## Train Dense Neural Network
using JLD2,Parameters, Flux, Random
using Flux:params, onehotbatch, crossentropy, onecold, throttle

N = length(Y)
permutation = randperm(N)
X_shuffled = X[:,permutation]
Y_shuffled = Y[permutation]
Nepochs = 1000
losses = zeros(Nepochs)

Ntrain = round(Int,0.8*N)
    Xtrain = X_shuffled[:,1:Ntrain]
    Ytrain = onehotbatch(Y_shuffled[1:Ntrain], mus)

    NN = Chain(
    Dense(L*L,100, relu),
    Dense(100, 100),
    Dense(100, length(mus)),
    softmax)

    opt = Adam()
    loss(X, y) = crossentropy(NN(X), y)
    progress = () -> @show(loss(X, y)) # callback to show loss

    z = @elapsed for i in 1:Nepochs
        println("ep = $i/$Nepochs")
        Flux.train!(loss, Flux.params(NN),[(Xtrain,Ytrain)], opt)#,cb = throttle(progress, 10))
        losses[i] = loss(Xtrain, Ytrain)
    end

prinz(z)
# Check whether the NN has at least learnt the train set
resultats = [abs(onecold(NN(Xtrain[:,i])) - onecold(Ytrain[:,i])) ≤ 2 for i in 1:Ntrain]
    mean(resultats)

resultats = [arclength(mus[onecold(NN(Xtrain[:,i]))],mus[onecold(Ytrain[:,i])],2pi) for i in 1:Ntrain]
    mean(resultats)

plot(losses,axis=:log,m=:circle)
    plot!(x->30x^(-0.5))
    # plot!(x->8x^(-0.2))

histogram([onecold(NN(Xtrain[:,i])) - onecold(Ytrain[:,i]) for i in 1:Ntrain],normalize=true)

# Visualize it
model = XY(params)
lattice = SquareLattice(L)
ind = rand(1:Ntrain)
ind = rand(findall(x->abs(x)>0.3,resultats))
    prediction = (onecold(NN((Xtrain[:,ind]))) == onecold(Ytrain[:,ind]))
    thetass = reshape(Xtrain[:,ind],2WINDOW + 1,2WINDOW + 1)
    titre = string(mus[onecold(Ytrain[:,ind])])*" vs "*string(mus[onecold(NN((Xtrain[:,ind])))])*" truth vs pred"
    p = plot_thetas(thetass,model,lattice,title=titre)
    display_quiver!(p,thetass,WINDOW)

# Testset
Ntest = N - Ntrain
Xtest = X_shuffled[:,Ntrain+1:end]
Ytest = onehotbatch(Y_shuffled[Ntrain+1:end], mus)

resultats = [abs(onecold(NN(Xtest[:,i])) - onecold(Ytest[:,i])) ≤ 2 for i in 1:Ntest]
    mean(resultats)

resultats = [arclength(mus[onecold(NN(Xtest[:,i]))],mus[onecold(Ytest[:,i])],2pi) for i in 1:Ntest]
    mean(resultats)

histogram([onecold(NN(Xtest[:,i])) - onecold(Ytest[:,i]) for i in 1:Ntest],normalize=true)


## Creation of actual systems for testing the NN efficiency on in vivo defects
include(srcdir("../parameters.jl"));
model = SPP(params) ; lattice = TriangularLattice(L)
thetas = init_thetas(model,lattice,params_init=params_init)
tmax = 2000
update!(thetas,model,lattice,tmax=tmax)
    plot_thetas(thetas,model,lattice,defects=false)
delta = 100
update!(thetas,model,lattice,delta)
    plot_thetas(thetas,model,lattice,defects=false)
number_defects(thetas,model,lattice)

zoom_quiver(thetas,model,lattice,20,95)
plot_thetas(thetas,model,lattice,defects=true)

# jldsave("data/for_ML/thetas/L$(L)_$(symmetry)Moving_rho$(rho)_A$(A)_tmax$(tmax)_n$(number_defects(thetas,model,lattice)).jld2";thetas)
def = spot_defects(thetas,model,lattice)
ind += 1
    i,j = def[1][ind][1:2]
    zoom_quiver(thetas,model,lattice,i,j)



jldsave
## Test in vivo
load("data/for_ML/")


## Augmentation of the base_dataset for CNN
using Augmentor
model = XY(params) ; lattice = TriangularLattice(L,periodic=false)
rot_angles = 0:15:270
    nb_noise = 200
    N = length(rot_angles)*nb_noise*length(mus)
    X = zeros(Float32,L,L,1,N)
    Y = zeros(Float32,N)
    token = 1
    for i in each(mus)
        mu = mus[i]
        for rotation_angle in rot_angles
            for n in 1:nb_noise
                ppl = Rotate(rotation_angle) |> Resize(L,L)
                augmented_thetas = augment(base_dataset[:,:,i],ppl) .+ Float32(deg2rad(rotation_angle)) + Float32(0.2)*rand(Float32)*randn(Float32,L,L)
                # relax!(augmented_thetas,model,0.5)
                X[:,:,1,token] = augmented_thetas
                Y[token] = mu
                token += 1
            end
        end
    end
X
ind = rand(1:size(X,4))
    p=plot_thetas(X[:,:,1,ind],model,lattice,defects=false)
    display_quiver!(p,X[:,:,1,ind],WINDOW)

Y
## Train Convolutional Neural Network
using JLD2,Parameters, Flux, Random
using Flux:params, onehotbatch, crossentropy, onecold, throttle

N = length(Y)
permutation = randperm(N)
X_shuffled = X[:,:,:,permutation]
Y_shuffled = Y[permutation]
Nepochs = 10
losses = zeros(Nepochs)

Ntrain = round(Int,0.8*N)
    Xtrain = X_shuffled[:,:,:,1:Ntrain]
    Ytrain = onehotbatch(Y_shuffled[1:Ntrain], mus)

    NN = Chain(
        Conv((3, 3), 1=>16, pad=(1,1), relu),
        x -> maxpool(x, (2,2)),
        Conv((3, 3), 16=>32, pad=(1,1), relu),
        x -> maxpool(x, (2,2)),
        # Reshape 3d tensor into a 2d one, at this point it should be (3, 3, 32, N)
        # hence the 288 (=3x3x32) in the `Dense` layer below:
        x -> reshape(x, :, size(x, 4)),
        Dense(288, length(mus)),
        softmax)


    opt = Adam()
    loss(X, y) = crossentropy(NN(X), y)
    progress = () -> @show(loss(X, y)) # callback to show loss

    z = @elapsed for i in 1:Nepochs
        println("ep = $i/$Nepochs")
        Flux.train!(loss, Flux.params(NN),[(Xtrain,Ytrain)], opt)#,cb = throttle(progress, 10))
        losses[i] = loss(Xtrain, Ytrain)
    end
    prinz(z)

prinz(100z)
# Check whether the NN has at least learnt the train set
resultats = abs.(onecold(NN(Xtrain[:,:,:,:])) - onecold(Ytrain)) .≤ 2
    mean(resultats)


resultats = arclength.(mus[onecold(NN(Xtrain))],mus[onecold(Ytrain)],2pi)
    mean(resultats)

histogram(onecold(NN(Xtrain)) - onecold(Ytrain),normalize=true)

plot(losses,axis=:log)
    plot!(x->8x^(-0.2))
    plot!(x->30x^(-0.45))

# Visualize it
ind = rand(1:Ntrain)
# ind = rand(findall(x->x==false,resultats))
    prediction = (onecold(NN((Xtrain[:,ind]))) == onecold(Ytrain[:,ind]))
    thetass = reshape(Xtrain[:,ind],2WINDOW + 1,2WINDOW + 1)
    titre = string(mus[onecold(Ytrain[:,ind])])*" vs "*string(mus[onecold(NN((Xtrain[:,ind])))])
    p = plot_thetas(thetass,model,lattice,title=titre)
    display_quiver!(p,thetass,WINDOW)


# Testset
Xtest = X_shuffled[:,:,:,Ntrain+1:end]
Ytest = onehotbatch(Y_shuffled[Ntrain+1:end], mus)

resultats = abs.(onecold(NN(Xtest)) - onecold(Ytest)) .≤ 3
    mean(resultats)
histogram(onecold(NN(Xtest)) - onecold(Ytest),normalize=true)
mus[onecold(NN(reshape(rand(15,15),(15,15,1,1))))...]


## Test in vivo and compare to

## Generate single defects and pairs of defects to verify coherence
include(srcdir("../parameters.jl"));

model = XY(params)
lattice = SquareLattice(L)
params_init["type1defect"] = "32"
params_init["type2defect"] = ["source","random"]
params_init["type2defect"] = "random"

    thetas = init_thetas(model,lattice,params_init=params_init)
    p=plot_thetas(thetas,model,lattice,defects=false)
    display_quiver!(p,thetas,9)

## Some variety of +1 defects
include(srcdir("../parameters.jl"));
model = XY(params)
lattice = SquareLattice(L)

mus = collect(Float16,0:0.4:2pi)
plotss = Vector(undef,length(mus))
for i in each(mus)
    params_init["type1defect"] = mus[i]
    thetas = init_thetas(model,lattice,params_init=params_init)
    p=plot_thetas(thetas,model,lattice,defects=false,colorbar=false,title="µ=$(mus[i])")
    display_quiver!(p,thetas,7)
    plotss[i] = p
end
plot(plotss...,layout=(4,4),size=(1600,1600))
# savefig(plotsdir("illustration_defects/defects_1.png"))

## Create a dataset for later defects ML indentification
#= The goal here is, for each defect type (8 or 16, depending
on whether we take the +/- 1 defects), to create N noisy
configurations (T is random, below ~0.3). Around 85% of them
will be used for training, the rest for testing.

Each config should be a zoom around the defect, a 11x11 square
centered on the defect.
=#
N = 1000 # the number of config for each defect
window = 7 # 7 x 7 square around the defect
# possible_defects = [(1/2,"source"),(1/2,"sink"),(1/2,"clockwise"),(1/2,"counterclockwise")]
# possible_defects = [(-1/2,"join"),(-1/2,"split"),(-1/2,"threefold1"),(-1/2,"threefold2")]
# possible_defects = [(1,"source"),(1,"sink")]
possible_defects = [(1,"source"),(1,"sink"),(1,"clockwise"),(1,"counterclockwise")]
# possible_defects = [(-1,"join"),(-1,"split"),(-1,"threefold1"),(-1,"threefold2")]
# possible_defects = [(1,"source"),(1,"sink")]
params["init"] = "single"
params["symmetry"] = "polar"
params["L"] = 80
X = zeros(Float32,2window+1,2window+1,N*length(possible_defects))
Y = vcat([fill(possible_defects[i][2],N) for i in 1:length(possible_defects)]...)

z = @elapsed for k in 1:length(possible_defects)
    println("Simulation for ",possible_defects[k]," started.")
    defect = possible_defects[k]
    params["q"] , params["type1defect"] = defect
    lattice = TriangularLattice(params["L"],periodic=false)
    model = XY(params)
    params["q"] > 0 ? ind = 1 : ind = 2
    for n in 1:N
        model.T = 0.3*rand()
        model.t = 0
        thetas = randomly_rotate(init_thetas(lattice,params=params))
        update!(thetas,model,lattice,50)
        i,j = spot_defects(thetas,model,lattice,find_type=false)[ind][1][1:2]
        no_problem,thetas_zoom = zoom(thetas,lattice,i,j,window)
        if no_problem  X[:,:,(k-1)*N + n] = thetas_zoom
        else X[:,:,(k-1)*N + n] = randomly_rotate(init_thetas(lattice,params=params))
        end
    end
end
prinz(z)
ind = rand(1:size(X,3))
    p = plot_thetas(X[:,:,ind],model,lattice,defects=false)
    display_quiver!(p,X[:,:,ind],window)
    xlims!(1,2window+1) ; ylims!(1,2window+1)

using JLD2
# jldsave(datadir("for_ML/dataset_T0.3Random_negative12defects_N$(N)_W$(window).jld2");X,Y,N,window,possible_defects,params,comments="Evolution time = 4, dt = 1E-2, Float32.")
# jldsave(datadir("for_ML/dataset_T0.3Random_negative1defects_N$(N)_W$(window).jld2");X,Y,N,window,possible_defects,params,comments="Evolution time = 4, dt = 1E-2, Float32.")

## See whether I can learn the features from a Dense Neural Network
using JLD2,Parameters, Flux, Random
using Flux:params, onehotbatch, crossentropy, onecold, throttle

# @unpack X,Y,window,possible_defects,comments,N = load(datadir("for_ML/dataset_T0.3Random_all12defects_N2000_W9.jld2"))
@unpack X,Y,window,possible_defects,comments,N = load(datadir("for_ML/dataset_T0.3Random_positive1defects_N1000_W7.jld2"))
permutation = randperm(size(X,3))
X_shuffled = X[:,:,permutation]
Y_shuffled = Y[permutation]
possible_labels = unique(Y)
Nepochs = 500
NN = 0
Ntrain = round(Int,0.8*length(Y))
    L = 2window + 1
    Xtrain = zeros(L*L,Ntrain) ; for i in 1:Ntrain  Xtrain[:,i] = vec(X_shuffled[:,:,i]) end
    Ytrain = onehotbatch(Y_shuffled[1:Ntrain], possible_labels)

    NN = Chain(
    Dense(L*L,40, relu),
    Dense(40, length(possible_labels)),
    softmax)

    opt = Adam()
    loss(X, y) = crossentropy(NN(X), y)
    progress = () -> @show(loss(X, y)) # callback to show loss

    for i in 1:Nepochs
        Flux.train!(loss, Flux.params(NN),[(Xtrain,Ytrain)], opt)
    end


# Check whether the NN has at least learnt the train set
resultats = [onecold(NN(Xtrain[:,i])) == onecold(Ytrain[:,i]) for i in 1:Ntrain]
    mean(resultats)

# Visualize it
ind = rand(1:Ntrain)
# ind = rand(findall(x->x==false,resultats))
    prediction = (onecold(NN((Xtrain[:,ind]))) == onecold(Ytrain[:,ind]))
    thetass = reshape(Xtrain[:,ind],2window + 1,2window + 1)
    p = plot_thetas(thetass,model,lattice,title=possible_labels[onecold(Ytrain[:,ind])])
    # p = plot_thetas(thetass,model,lattice,title=possible_labels[onecold(Ytrain[:,ind])]*" vs "*possible_defects[onecold(NN((Xtrain[:,ind])))][2])
    display_quiver!(p,thetass,window)

## TestSet
Ntest  = length(Y) - Ntrain
    Xtest = zeros(L*L,Ntest)
    for i in 1:Ntest  Xtest[:,i] = vec(X_shuffled[:,:,Ntrain+i]) end
    Ytest = onehotbatch(Y_shuffled[Ntrain+1:end], possible_labels)

    resultats = [onecold(NN(Xtest[:,i])) == onecold(Ytest[:,i]) for i in 1:Ntest]
    mean(resultats)

# Visualize it
# model = XY(params)
# lattice = TriangularLattice(L,periodic=false)
ind = rand(1:Ntest)
ind = rand(findall(x->x==false,resultats))
    prediction = (onecold(NN((Xtest[:,ind]))) == onecold(Ytest[:,ind]))
    thetass = reshape(Xtest[:,ind],2window + 1,2window + 1)
    p = plot_thetas(thetass,model,lattice,title=possible_labels[onecold(Ytest[:,ind])]*" vs "*possible_defects[onecold(NN((Xtest[:,ind])))][2])
    display_quiver!(p,thetass,window)

## Save NN
# cd("D:/Documents/Research/projects/LatticeModels")
# comments="Evolution time = 4, dt = 1E-2, Float32, T<0.3 random, $(N) config per defect, originally from 32x32 lattices, $(Nepochs) epochs to train NN "
# jldsave("NN_negative_1_defects_N$(N)_W$(window).jld2";NN,possible_defects=possible_labels,comments=comments)

## Test it on completely new images
using Flux
using Flux:onecold
NN = load(datadir("for_ML/NN_all_12_defects.jl"),"NN")
NN(rand(121))
onecold(NN(rand(121)))

# Generate new defect
lattice = TriangularLattice(L)
model = XY(params)
thetas = init_thetas(lattice,params=params)
# update!(thetas,model,lattice,200)
p=plot_thetas(thetas,model,lattice) # plot the whole field, to then zoom on specific defects

i,j = (118,35) # loc where zoom
window = 5
relax!(thetas,model,0.1)
thetas_zoom = (thetas[i-window:i+window,j-window:j+window])
    p=plot_thetas(thetas_zoom,model,lattice)
    display_quiver!(p,thetas_zoom,window)
onecold(NN(vec(thetas_zoom))) # check whether its corresponds to reality
# @btime onecold(NN(vec(thetas_zoom))) # 1.5 µs (including 0.05 µs for onecold(...))

# update for a little time at "high temperature" T = 0.4 to see whether the algo is robust
model.T = 0.4
update!(thetas,model,lattice,202)

# Conclusion : le champ theta est très brouillon (even at T small for rho = 1) mais ca a l'air de fonctionner

## Different -1 defects (with L = 32, single)
# Result, already known: split<->join and threefold1<->threefold2 are related via a 90° rotation. So the NN cannot really learn the difference as there is no.
# include(srcdir("../parameters.jl"));
# lattice = TriangularLattice(L,periodic=false)
# model = XY(params)
#
# params["type1defect"] = "split"
# thetas = init_thetas(lattice,params=params)
# zoom_quiver(thetas,model,lattice,17,17,10) ; title!("split")
#
# params["type1defect"] = "join"
# thetas = init_thetas(lattice,params=params)
# zoom_quiver(rotate_clockwise90(thetas),model,lattice,17,17,10) ; title!("rotated join = split")
#
# params["type1defect"] = "threefold1"
# thetas = init_thetas(lattice,params=params)
# zoom_quiver(rotate_180(thetas),model,lattice,17,17,10) ; title!("rotated ccw = split")

## Check the training set, maybe the errors in classification comes from badly labelled data in the beggining


## Plot the different defects and defect pairs
include(srcdir("../parameters.jl"));
    cols = cgrad([:black,:blue,:green,:orange,:red,:black])
    params["L"] = 20
    window = Int(params["L"]/2-1)
    model = XY(params) # in fact, useless for plotting defects at t = 0
    lattice = SquareLattice(L,periodic=true,single=true) # in fact, useless for plotting defects at t = 0

# +1 Defects
plotsP1 = []
    params["q"] = +1
    params["symmetry"] = "polar"
    params["init"] = "single"
    types = ["source","sink","clockwise","counterclockwise"]
    for type in types
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="+1 "*type)
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsP1,p)
    end
    pP1 = plot(plotsP1...,layout=(1,4),size=(400*4,400))

# -1 Defects
plotsM1 = []
    params["q"] = -1
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["convergent","divergent","threefold1","threefold2"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="-1",size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsM1,p)
    end
    pM1 = plot(plotsM1...,layout=(1,4),size=(400*4,400))

# +1/2 Defects
plotsP12 = []
    params["q"] = +1/2
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["source","sink","clockwise","counterclockwise"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="+1/2 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsP12,p)
    end
    pP12 = plot(plotsP12...,layout=(1,4),size=(400*4,400))

# -1/2 Defects
plotsM12 = []
    params["q"] = -1/2
    params["symmetry"] = "polar"
    params["init"] = "single"
    for type in ["join","split","threefold1","threefold2"]
        params["type1defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="-1/2 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plotsM12,p)
    end
    pM12 = plot(plotsM12...,layout=(1,4),size=(400*4,400))

# 4 ≠ Pairs (all the others are equivalent to one of those)
plots_pairs = []
    params["r0"] = 8
    params["q"] = 1/2
    params["symmetry"] = "polar"
    params["init"] = "pair"
    for type in ["pair1","pair2","pair3","pair4"]
        params["type2defect"] = type
        thetas = init_thetas(lattice,params=params)
        p = plot_thetas(thetas,model,lattice,colorbar=false,title="1/2 "*type,size=(400,400))
        display_quiver!(p,thetas,window)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
        push!(plots_pairs,p)
    end
    pPairs = plot(plots_pairs...,layout=(1,4),size=(1600,400))

# savefig(pP1,plotsdir("illustration_defects/P1_defects.png"))
# savefig(pP1,plotsdir("illustration_defects/P1_defects.svg"))
#
# savefig(pP12,plotsdir("illustration_defects/P12_defects.png"))
# savefig(pP12,plotsdir("illustration_defects/P12_defects.svg"))
#
# savefig(pM1,plotsdir("illustration_defects/M1_defects.png"))
# savefig(pM1,plotsdir("illustration_defects/M1_defects.svg"))
#
# savefig(pM12,plotsdir("illustration_defects/M12_defects.png"))
# savefig(pM12,plotsdir("illustration_defects/M12_defects.svg"))
#
# savefig(pPairs,plotsdir("illustration_defects/Pairs_defects.png"))
# savefig(pPairs,plotsdir("illustration_defects/Pairs_defects.svg"))

## Discriminate between same q but different defects : divergence and rotationnal ?
# First implement, then test on my defects, then export to recognition by DefectTracker
include(srcdir("../parameters.jl"));
params["type1defect"] = "source"
    params["q"] = 1/2
    lattice = SquareLattice(L)
    thetas = init_thetas(lattice,params=params)
    model = XY(params)
    window = 9 # for L = 20
    divergence,rotational = get_div_rot(thetas,lattice)

p = plot_thetas(thetas,model,lattice,defects=true)
    display_quiver!(p,thetas,window)
    xlims!(1,2window+1) ; ylims!(1,2window+1)
    pos = spot_defects(thetas,model,lattice)[1][1][1:2]
    scatter!(pos,m=:xcross,c=:white)

heatmap(divergence',aspect_ratio=1,size=(485,450),c=cgrad([:blue,:white,:red]),title="div")
    pos = spot_defects(thetas,model,lattice)[1][1][1:2]
    scatter!(pos,m=:xcross,c=:black)
heatmap(rotational',aspect_ratio=1,size=(485,450),c=cgrad([:blue,:white,:red]),title="rot")
    pos = spot_defects(thetas,model,lattice)[1][1][1:2]
    scatter!(pos,m=:xcross,c=:black)

get_divergence(get_neighbours(thetas,model,lattice,10,10))
get_rotational(get_neighbours(thetas,model,lattice,10,10))

## Now observe them in a noisy environnement
include(srcdir("../parameters.jl"));
params["type1defect"] = "join"
    params["q"] = -1/2
    params["T"] = 0.2
    lattice = TriangularLattice(L,)
    model = XY(params)
    thetas = init_thetas(lattice,params=params)
    update!(thetas,model,lattice,5)
    # relax!(thetas,model,0.5) # because in real use, theta will be relaxed for t=0.3
    window = 9 # for L = 20
    divergence,rotational = get_div_rot(thetas,lattice)

    params["q"] > 0 ? ind = 1 : ind = 2
    p = plot_thetas(thetas,model,lattice,defects=false)
    display_quiver!(p,thetas,window)
    xlims!(1,2window+1) ; ylims!(1,2window+1)
pos = spot_defects(thetas,model,lattice)[ind][1][1:2]
    scatter!(pos,m=:xcross,ms=7,c=:white)
    pd = heatmap(divergence',aspect_ratio=1,size=(485,400),c=cgrad([:blue,:white,:red]))
        pos = spot_defects(thetas,model,lattice)[ind][1][1:2]
        scatter!(pos,m=:xcross,ms=7,c=:black)
        xlims!(1,2window+1) ; ylims!(1,2window+1)
    pr = heatmap(rotational',aspect_ratio=1,size=(485,400),c=cgrad([:blue,:white,:red]))
    pos = spot_defects(thetas,model,lattice)[ind][1][1:2]
    scatter!(pos,m=:xcross,ms=7,c=:black)
    xlims!(1,2window+1) ; ylims!(1,2window+1)
    plot(p,pd,pr,layout=(3,1),size=(485,1200))
