using Parameters

# export XY, ForcedXY, VisionXY

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
    rho::AbstractFloat

end
function XY(params)
    @unpack T,symmetry,dt,float_type,rho = params
    T,dt,rho = convert.(float_type,(T,dt,rho))

    return XY{float_type}(T,symmetry,dt,zero(float_type),rho)
end


## ------------------------- Forced / Active XY Model -------------------------
mutable struct ForcedXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    const omegas::Matrix{AbstractFloat}
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function ForcedXY(params)
    @unpack L,T,Var,symmetry,dt,float_type,rho  = params
    T,Var,dt,rho = convert.(float_type,(T,Var,dt,rho))

    omegas = sqrt(Var)*randn(float_type,L,L)

    return ForcedXY{float_type}(T,Var,symmetry,omegas,dt,zero(float_type),rho)
end

## ------------------------- Moving XY -------------------------
mutable struct MovingXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    A::AbstractFloat
    symmetry::String
    propulsion::String
    t::AbstractFloat
    rho::AbstractFloat
    algo::String
    width_proposal::AbstractFloat # to be hidden from user once benchmarked
end
function MovingXY(params)
    @unpack T,A,rho,symmetry,algo,propulsion,float_type,width_proposal = params
    T,A,rho,width_proposal = convert.(float_type,(T,A,rho,width_proposal))

    return MovingXY{float_type}(T,A,symmetry,propulsion,zero(float_type),rho,algo,width_proposal)
end


## ------------------- Non Reciprocal (Vision Cone) XY Model -------------------
mutable struct VisionXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    vision::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function VisionXY(params)
    @unpack T,vision,symmetry,dt,float_type,rho = params
    T,vision,dt,rho = convert.(float_type,(T,vision,dt,rho))

    if vision ≠ 2π @assert symmetry == "polar" "I am not sure how to interpret a vision cone with nematic symmetry" end

    return VisionXY{float_type}(T,vision,symmetry,dt,zero(float_type),rho)
end
