using DrWatson ; @quickactivate "LatticeModels"
include(srcdir("LatticeModels.jl"))
using Plots,ColorSchemes,LaTeXStrings
pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()

include(srcdir("../parameters.jl"));

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
# possible_defects = [(1,"source"),(1,"sink"),(1,"clockwise"),(1,"counterclockwise")]
possible_defects = [(-1,"join"),(-1,"split"),(-1,"threefold1"),(-1,"threefold2")]
# possible_defects = [(1,"source"),(1,"sink")]
params["init"] = "single"
params["symmetry"] = "polar"
params["L"] = 32
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
        update!(thetas,model,lattice,4)
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
# permutation = randperm(size(X,3))
# X_shuffled = X[:,:,permutation]
# Y_shuffled = Y[permutation]
# possible_labels = unique(Y)
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
    prediction = (onecold(NN(Xtrain[:,ind])) == onecold(Ytrain[:,ind]))
    thetass = reshape(Xtrain[:,ind],2window + 1,2window + 1)
    p = plot_thetas(thetass,model,lattice,title=possible_labels[onecold(Ytrain[:,ind])]*" , "*string(prediction))
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

## Different -1 defects
include(srcdir("../parameters.jl"));
lattice = TriangularLattice(L,periodic=false)
model = XY(params)
thetas = init_thetas(lattice,params=params)
zoom_quiver(thetas,model,lattice,17,17,10) ; title!("split")

include(srcdir("../parameters.jl"));
lattice = TriangularLattice(L,periodic=false)
model = XY(params)
thetas = init_thetas(lattice,params=params)
zoom_quiver(rotate_clockwise90(thetas),model,lattice,17,17,10) ; title!("split")

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
