logspace(x1, x2, n) = [10.0 ^y for y in range(log10(x1), log10(x2), length=n)]
prinz(z) = println("Runtime : $(round(Int,z)) seconds = $(round(z/60,digits=2)) minutes = $(round(z/3600,digits=2)) hours.")
each = eachindex

nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,y) = mapslices(nanstd,x,dims=y)
replace_nan_with_zeros(v) = map(x -> isnan(x) ? zero(x) : x, v)


import Base: +,-
+(x::Vector{Tuple{T,T}},y::Vector{Tuple{T,T}}) where T<:Number = [x[i] .+ y[i] for i in 1:length(x)]
-(x::Vector{Tuple{T,T}},y::Vector{Tuple{T,T}}) where T<:Number = [x[i] .- y[i] for i in 1:length(x)]
+(x::Vector{Tuple{T,T}},y::Tuple{T,T}) where T<:Number = [x[i] .+ y for i in 1:length(x)]
+(y::Tuple{T,T},x::Vector{Tuple{T,T}}) where T<:Number = [x[i] .+ y for i in 1:length(x)]
-(x::Vector{Tuple{T,T}},y::Tuple{T,T}) where T<:Number = [x[i] .- y for i in 1:length(x)]
-(y::Tuple{T,T},x::Vector{Tuple{T,T}}) where T<:Number = [x[i] .- y for i in 1:length(x)]

function index_element_repeated(x)
    indices = []
    L = length(x)
    for i in 1:L
        count = 0
        for j in 1:L
            if x[j] == x[i] count += 1 end
        end
        if count > 1 push!(indices,i) end
    end
    return indices
end

function remove_negative(input)
    array = Float32.(copy(input))
    for i in 1:length(array)
        if array[i] ≤ 0 array[i] = NaN end
    end
    return array
end

function corr_length(rs,C::Vector{T},seuil=exp(-1))::T where T<:AbstractFloat # from a time series, returns the correlation length ()
    i_after = findfirst(x->x<seuil,C)
    if i_after ≠ nothing && i_after > 1
        # Linear interpolation
        i_befor = i_after - 1
        r_after = rs[i_after]
        r_befor = rs[i_befor]
        c_after = C[i_after]
        c_befor = C[i_befor]
        ξ = (seuil*(r_after-r_befor) -  (c_befor*r_after - r_befor*c_after))/(c_after-c_befor)
    else
    ξ = NaN
    end
    return ξ
end

function exp_moving_average(x,window)
    y = zeros(length(x))
    for t in 1:length(y)
        y[t] = sum([x[t-i]*((window-1)/window)^(i) for i in 0:t-1])/window
    end
    return y
end

function smooth(X;over=3) ## for smoother plots
    smoothed = copy(X)
    coeff = [2^(i-1) for i in 1:over]
    coeff = coeff./sum(coeff)
    s = length(coeff)
    for i in 1+s:length(smoothed)
        smoothed[i] = X[i-s+1:i]'*coeff
    end
    return smoothed
end
