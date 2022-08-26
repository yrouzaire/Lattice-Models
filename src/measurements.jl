include("lattice_general.jl");
include("models.jl");

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
