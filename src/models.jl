export XY, AXY, VisionXY

abstract type AbstractModel{AbstractFloat} end
abstract type AbstractModel2 end
function sym(model::AbstractModel{T}) where T<:AbstractFloat
    model.symmetry == "polar" ? symm = 2π : symm = π
    return T(symm)
end

## ---------------------------- Classical XY Model ----------------------------
mutable struct XY{AbstractFloat} <: AbstractModel{AbstractFloat}
    L::Int
    T::AbstractFloat
    symmetry::String
    thetas::Matrix{AbstractFloat}
end
function XY(L::Int,T,symmetry::String;float_type=Float32)
    thetas = zeros(float_type,L,L)
    T = convert(float_type,T)

    return XY{float_type}(L,T,symmetry,thetas)
end


## ------------------------- Forced / Active XY Model -------------------------
mutable struct AXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    L::Int
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    thetas::Matrix{AbstractFloat}
    const omegas::Matrix{AbstractFloat}
end
function AXY(L::Int,T,Var,symmetry::String;float_type=Float32)
    thetas = zeros(float_type,L,L)
    T,Var = convert.(float_type,(T,Var))
    omegas = sqrt(Var)*randn(L,L)

    return AXY{float_type}(L,T,Var,symmetry,thetas,omegas)
    # if Var == 0
    #     return XY{float_type}(L,T,symmetry,thetas)
    # else
    #     return AXY{float_type}(L,T,Var,symmetry,thetas,omegas)
    # end
end

## ------------------- Non Reciprocal (Vision Cone) XY Model -------------------
mutable struct VisionXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    L::Int
    T::AbstractFloat
    vision::AbstractFloat
    symmetry::String
    thetas::Matrix{AbstractFloat}
end
function VisionXY(L::Int,T,vision,symmetry::String;float_type=Float32)
    thetas = zeros(float_type,L,L)
    T,vision = convert.(float_type,(T,vision))

    return VisionXY{float_type}(L,T,vision,symmetry,thetas)
    # if vision == 2π
    #     return XY{float_type}(L,T,symmetry,thetas)
    # else
    #     return VisionXY{float_type}(L,T,vision,symmetry,thetas)
    # end
end
