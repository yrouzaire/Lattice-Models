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

function infer_mu(thetas;q,window=WINDOW,decay=true)
    if decay && q > 0 return infer_mu_decay(thetas,q=q,window=window)
    else return infer_mu_0(thetas,q=q,window=window)
    end
end
function infer_mu_0(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    # tmp = Complex(0)
    range = 2:L-1
    for j in range, i in range
        muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
        # i<->j irrelevant because i,j and j,i have the same weight for "mean" operation
    end
    moyenne = angle(mean(exp.(im*muss[range,range])))
    if     q == 1    correction = pi - 0.33228605
    elseif q == -1   correction = pi/2 + 0.04
    elseif q == 1/2  correction = 0.1
    elseif q == -1/2 correction = pi/4
    end
    #= To be honest, I don't know where the shifts might come from,
    In the beggining, I thought maybe from the spin at the center of the defect,
    where theta should not be defined. But if one changes mean->sum and adds the condition
    "if (i == j == window)", the result only becomes weirder... =#
    # return mod.(muss,2pi)
    return mod(moyenne .+ correction,2π)
end
function infer_mu_decay(thetas::Matrix{T};q,window=WINDOW) where T<:AbstractFloat
    L = size(thetas,1)
    @assert L == 2window+1
    muss = zeros(size(thetas))
    tmp = Complex(0)
    range = 2:L-1
    for j in range, i in range
        muss[i,j] = thetas[i,j] - abs(q)*atan( (i-window) ,(j-window))
        # i<->j irrelevant because i,j and j,i have the same weight for "mean" or "sum" operation
        distance = sqrt((i-window)^2 + (j-window)^2)
        # if distance == 0 distance = Inf end # does horribly wrong
        tmp += exp(im*muss[i,j] - 1*distance) # 1*distance seems to be the best
    end
    moyenne = angle(tmp)
    if     q == 1   correction = pi - 0.603228605 # 0.33228605 was the original cst
    elseif q == -1  correction = pi/2 + 0.18
    elseif q == 1/2 correction = 0.1
    elseif q == -1/2 correction = pi/4
    end
    #= To be honest, I don't know where the shifts might come from,
    In the beggining, I thought maybe from the spin at the center of the defect,
    where theta should not be defined. But if one changes mean->sum and adds the condition
    "if (i == j == window)", the result only becomes weirder... =#
    return mod(moyenne .+ correction,2π)
end

## Test noiseless
# inferredP  = zeros(length(mus))
#     inferredN = zeros(length(mus))
#     decay = true
#     for ind in 1:64
#         inferredP[ind] = infer_mu(base_datasetP[:,:,ind],q=CHARGE,decay=decay)
#         inferredN[ind] = infer_mu(base_datasetN[:,:,ind],q=-CHARGE,decay=decay)
#     end
#     p1 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noiseless")
#     plot!(mus,inferredP,m=true,c=:red,label="q=1")
#     plot!(mus,inferredN,m=true,c=:green,label="q=-1")
#     plot!(x->x,c=:black)

## Test with noise
# inferredP  = zeros(length(mus))
#     inferredN = zeros(length(mus))
#     noise = 0.3randn(W21,W21,64)
#     decay = false
#     for ind in 1:64
#         inferredP[ind]  = infer_mu(base_datasetP[:,:,ind] + noise[:,:,ind],q=CHARGE,decay=decay)
#         inferredN[ind] = infer_mu(base_datasetN[:,:,ind]+ noise[:,:,ind],q=-CHARGE,decay=decay)
#     end
#     p2 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noisy (without decay)")
#     plot!(mus,inferredP,line=false,m=true,c=:red,label="q=1")
#     plot!(mus,inferredN,line=false,m=true,c=:green,label="q=-1")
#     plot!(x->x,c=:black)
#
# decay = true
#     for ind in 1:64
#         inferredP[ind]  = infer_mu(base_datasetP[:,:,ind] + noise[:,:,ind],q=CHARGE,decay=decay)
#         inferredN[ind] = infer_mu(base_datasetN[:,:,ind]+ noise[:,:,ind],q=-CHARGE,decay=decay)
#     end
#     p3 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noisy (with decay for q=1)")
#     plot!(mus,inferredP,line=false,m=true,c=:red,label="q=1")
#     plot!(mus,inferredN,line=false,m=true,c=:green,label="q=-1")
#     plot!(x->x,c=:black)
#
#
# plot(p1,p2,p3,size=(1200,400),layout=(1,3))
# # savefig("plots/procedure_infer_mu_polar.png")
