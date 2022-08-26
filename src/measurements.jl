include("lattice_general.jl");
include("models.jl");


function OP(thetas::Matrix{T})::Vector{T} where T<:AbstractFloat
    tmp_polar   = Complex(0)
    tmp_nematic = Complex(0)
    for theta in thetas
        if !isnan(theta)
            tmp_polar   += exp(im*theta)
            tmp_nematic += exp(2im*theta)
        end
    end
    return abs(tmp_polar/length(thetas)) , abs(tmp_nematic/length(thetas))
end

function local_mag(thetas::Matrix{T},L::Int,l::Int)::Matrix{T} where T<:AbstractFloat
    # L is the size of the system
    # l is the size of the local box for averaging the magnetisation
    @assert isinteger(L/l)
    c = Int(L/l)
    δ = zeros(2,c^2)
    for j in 1:c
        for i in 1:c
            δ[:,square_to_linear_index(i,j,c)] = get_Delta(thetas[1+l*(i-1):i*l, 1+l*(j-1):j*l])
        end
    end
    return δ
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
