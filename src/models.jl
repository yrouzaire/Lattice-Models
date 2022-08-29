using Parameters

export XY, AXY, VisionXY

abstract type AbstractModel{AbstractFloat} end
function sym(model::AbstractModel{T}) where T<:AbstractFloat
    model.symmetry == "polar" ? symm = 2π : symm = π
    return T(symm)
end

## ---------------------------- Classical XY Model ----------------------------
mutable struct XY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    dt::AbstractFloat
end
function XY(params_phys,params_num)
    @unpack T,symmetry  = params_phys
    @unpack dt,float_type = params_num

    T,dt = convert.(float_type,(T,dt))

    return XY{float_type}(T,symmetry,dt)
end


## ------------------------- Forced / Active XY Model -------------------------
mutable struct AXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    const omegas::Matrix{AbstractFloat}
    dt::AbstractFloat
end
function AXY(params_phys,params_num)
    @unpack L,T,Var,symmetry  = params_phys
    @unpack dt,float_type = params_num
    T,Var,dt = convert.(float_type,(T,Var,dt))

    omegas = sqrt(Var)*randn(float_type,L,L)

    return AXY{float_type}(T,Var,symmetry,omegas,dt)
end


## ------------------- Non Reciprocal (Vision Cone) XY Model -------------------
mutable struct VisionXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    vision::AbstractFloat
    symmetry::String
    dt::AbstractFloat
end
function VisionXY(params_phys,params_num)
    @unpack T,vision,symmetry  = params_phys
    @unpack dt,float_type = params_num

    T,vision,dt = convert.(float_type,(T,vision,dt))

    if vision ≠ 2π @assert symmetry == "polar" "I am not sure how to interpret a vision cone with nematic symmetry" end

    return VisionXY{float_type}(T,vision,symmetry,dt)
end
