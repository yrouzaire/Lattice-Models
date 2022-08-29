using Parameters

export XY, AXY, VisionXY

abstract type AbstractModel{AbstractFloat} end
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
    thetas_new::Matrix{AbstractFloat}
    dt::AbstractFloat
end
function XY(params_phys,params_num)
    @unpack L,T,symmetry  = params_phys
    @unpack dt,float_type = params_num

    thetas = zeros(float_type,L,L)
    T,dt = convert.(float_type,(T,dt))

    return XY{float_type}(L,T,symmetry,thetas,thetas,dt) # thetas is written twice because thetas_new = thetas at t=0
end


## ------------------------- Forced / Active XY Model -------------------------
mutable struct AXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    L::Int
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    thetas::Matrix{AbstractFloat}
    thetas_new::Matrix{AbstractFloat}
    const omegas::Matrix{AbstractFloat}
    dt::AbstractFloat
end
function AXY(params_phys,params_num)
    @unpack L,T,Var,symmetry  = params_phys
    @unpack dt,float_type = params_num
    T,Var,dt = convert.(float_type,(T,Var,dt))

    thetas = zeros(float_type,L,L)
    omegas = sqrt(Var)*randn(float_type,L,L)

    return AXY{float_type}(L,T,Var,symmetry,thetas,thetas,omegas,dt) # thetas is written twice because thetas_new = thetas at t=0
end


## ------------------- Non Reciprocal (Vision Cone) XY Model -------------------
mutable struct VisionXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    L::Int
    T::AbstractFloat
    vision::AbstractFloat
    symmetry::String
    thetas::Matrix{AbstractFloat}
end
function VisionXY(params_phys,params_num)
    @unpack L,T,vision,symmetry  = params_phys
    @unpack float_type = params_num
    T,vision = convert.(float_type,(T,vision))

    if vision ≠ 2π @assert symmetry == "polar" "I am not sure how to interpret a vision cone with nematic symmetry" end

    thetas = zeros(float_type,L,L)
    return VisionXY{float_type}(L,T,vision,symmetry,thetas)
end
