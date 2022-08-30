using Parameters

export XY, ForcedXY, VisionXY

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
    t::AbstractFloat
end
function XY(params_phys,params_num)
    @unpack T,symmetry  = params_phys
    @unpack dt,float_type = params_num

    T,dt = convert.(float_type,(T,dt))

    return XY{float_type}(T,symmetry,dt,float_type(0))
end


## ------------------------- Forced / Active XY Model -------------------------
mutable struct ForcedXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    const omegas::Matrix{AbstractFloat}
    dt::AbstractFloat
    t::AbstractFloat
end
function ForcedXY(params_phys,params_num)
    @unpack L,T,Var,symmetry  = params_phys
    @unpack dt,float_type = params_num
    T,Var,dt = convert.(float_type,(T,Var,dt))

    omegas = sqrt(Var)*randn(float_type,L,L)

    return ForcedXY{float_type}(T,Var,symmetry,omegas,dt,float_type(0))
end

## ------------------------- Moving XY -------------------------
mutable struct MovingXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    A::AbstractFloat
    symmetry::String
    t::AbstractFloat
    rho::AbstractFloat
end
function MovingXY(params_phys,params_num)
    @unpack T,A,rho,symmetry  = params_phys
    @unpack float_type = params_num
    T,A,rho = convert.(float_type,(T,A,rho))

    return MovingXY{float_type}(T,A,symmetry,float_type(0),rho)
end


## ------------------- Non Reciprocal (Vision Cone) XY Model -------------------
mutable struct VisionXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    vision::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
end
function VisionXY(params_phys,params_num)
    @unpack T,vision,symmetry  = params_phys
    @unpack dt,float_type = params_num

    T,vision,dt = convert.(float_type,(T,vision,dt))

    if vision ≠ 2π @assert symmetry == "polar" "I am not sure how to interpret a vision cone with nematic symmetry" end

    return VisionXY{float_type}(T,vision,symmetry,dt,float_type(0))
end
