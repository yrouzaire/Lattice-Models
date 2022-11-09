cd("D:/Documents/Research/projects/LatticeModels")
 using DrWatson ; @quickactivate "LatticeModels"
 include(srcdir("LatticeModels.jl"))
 using Plots,ColorSchemes,LaTeXStrings
 pyplot(box=true,fontfamily="sans-serif",label=nothing,palette=ColorSchemes.tab10.colors[1:10],grid=false,markerstrokewidth=0,linewidth=1.3,size=(400,400),thickness_scaling = 1.5) ; plot()
&
@unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µP1.jld2")
    base_dataset1 = base_dataset
    @unpack base_dataset,mus,dµ = load("data/for_ML/base_dataset_µN12.jld2")
    base_dataset12 = base_dataset;

function infer_mu(thetas;q,window=WINDOW,decay=true)
    if decay return infer_mu_decay(thetas,q=q,window=window)
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
        # tmp += exp(im*muss[i,j] - 0.25sqrt((i-window)^2 + (j-window)^2))
    end
    # muss[window,window] = 0
    corrmat = zeros(size(thetas)) ; corrmat[window,window] = 1
    moyenne = angle(mean(exp.(im*muss[range,range])-0corrmat[range,range]))
    if     abs(q) == 1   correction = pi - 0.33228605 # 0.2 works perfectly with corrmat, 0.33228605 was the original cst
    elseif abs(q) == 1/2 correction = 0.8
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
        tmp += exp(im*muss[i,j] - sqrt((i-window)^2 + (j-window)^2))
    end
    moyenne = angle(tmp)
    if     q == 1   correction = pi - 0.603228605 # 0.33228605 was the original cst
    elseif q == 1/2 correction = 0.8
    elseif q == -1/2 correction = 0.25
    end
    #= To be honest, I don't know where the shifts might come from,
    In the beggining, I thought maybe from the spin at the center of the defect,
    where theta should not be defined. But if one changes mean->sum and adds the condition
    "if (i == j == window)", the result only becomes weirder... =#
    return mod(moyenne .+ correction,2π)
end

# ## Test noiseless
# inferred1  = zeros(length(mus))
# inferred12 = zeros(length(mus))
#     for ind in 1:64
#         inferred1[ind] = infer_mu(base_dataset1[:,:,ind],q=1,decay=false)
#         inferred12[ind] = infer_mu(base_dataset12[:,:,ind],q=-1/2,decay=false)
#     end
#     p1 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noiseless")
#     plot!(mus,inferred1,m=true,c=:red,label="Polar")
#     plot!(mus,inferred12,m=true,c=:green,label="Nematic")
#     plot!(x->x,c=:black)
#
# ## Test with noise
# inferred1  = zeros(length(mus))
# inferred12 = zeros(length(mus))
# noise = 0.3randn(W21,W21,64)
# decay = false
#     for ind in 1:64
#         inferred1[ind]  = infer_mu(base_dataset1[:,:,ind] + noise[:,:,ind],q=1,decay=decay)
#         inferred12[ind] = infer_mu(base_dataset12[:,:,ind]+ noise[:,:,ind],q=-1/2,decay=decay)
#     end
#     p2 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noisy (without decay)")
#     plot!(mus,inferred1,line=false,m=true,c=:red,label="Polar")
#     plot!(mus,inferred12,line=false,m=true,c=:green,label="Nematic")
#     plot!(x->x,c=:black)
#
# decay = true
#     for ind in 1:64
#         inferred1[ind]  = infer_mu(base_dataset1[:,:,ind] + noise[:,:,ind],q=1,decay=decay)
#         inferred12[ind] = infer_mu(base_dataset12[:,:,ind]+ noise[:,:,ind],q=-1/2,decay=decay)
#     end
#     p3 = plot(xlabel="True µ",ylabel="Inferred µ",legend=:top,title="Noisy (with decay)")
#     plot!(mus,inferred1,line=false,m=true,c=:red,label="Polar")
#     plot!(mus,inferred12,line=false,m=true,c=:green,label="Nematic")
#     plot!(x->x,c=:black)
#
#
# plot(p1,p2,p3,size=(1200,400),layout=(1,3))
# savefig("plots/procedure_infer_mu_polar.png")
