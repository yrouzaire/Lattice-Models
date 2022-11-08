using Parameters

abstract type AbstractModel{AbstractFloat} end
abstract type AbstractPropagationModel{AbstractFloat} <: AbstractModel{AbstractFloat} end

function modd(model::AbstractModel{FT}) where FT<:AbstractFloat
    if model.symmetry == "polar" return FT(2pi)
    elseif model.symmetry == "nematic" return FT(pi)
    end
end
function sym(model::AbstractModel{FT})::FT where FT<:AbstractFloat
    if model.symmetry == "polar" return 1.0
    elseif model.symmetry == "nematic" return 2.0
    end
end

## ---------------------------- Classical XY Model ----------------------------
abstract type AbstractXYModel{AbstractFloat} <: AbstractModel{AbstractFloat} end
function XY(params) # by default, return LangevinXY
    @unpack T,symmetry,dt,float_type,rho,algo = params
    if algo in ["MC","MonteCarlo"]
        return MonteCarloXY{float_type}(T,symmetry,zero(float_type),rho)
    elseif algo == "Langevin"
        return LangevinXY{float_type}(T,symmetry,dt,zero(float_type),rho)
    else
        println()
        println("WARNING ! Unknown algo provided for XY model: return XY with Langevin dynamics.")
        return LangevinXY{float_type}(T,symmetry,dt,zero(float_type),rho)
    end
end

mutable struct LangevinXY{AbstractFloat} <: AbstractXYModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat

end
function LangevinXY(params)
    @unpack T,symmetry,dt,float_type,rho = params
    T,dt,rho = convert.(float_type,(T,dt,rho))

    return LangevinXY{float_type}(T,symmetry,dt,zero(float_type),rho)
end

# MonteCarlo
mutable struct MonteCarloXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    symmetry::String
    t::AbstractFloat
    rho::AbstractFloat

end
function MonteCarloXY(params)
    @unpack T,symmetry,float_type,rho = params
    T,rho = convert.(float_type,(T,rho))

    return MonteCarloXY{float_type}(T,symmetry,zero(float_type),rho)
end


## ------------------------- Forced / Active XY Model -------------------------
mutable struct ForcedXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    omegas::Matrix{AbstractFloat}
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

## --------------------- Self Propelled Particles (SPP) ---------------------
mutable struct SPP{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    A::AbstractFloat
    symmetry::String
    propulsion::String
    t::AbstractFloat
    rho::AbstractFloat
    algo::String
end
function SPP(params)
    @unpack T,A,rho,symmetry,algo,propulsion,float_type = params
    T,A,rho = convert.(float_type,(T,A,rho))

    return SPP{float_type}(T,A,symmetry,propulsion,zero(float_type),rho,algo)
end


## ------------------- Non Reciprocal XY Models -------------------
# (Vision Cone)
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

# (Softly Tuned Couplings)
mutable struct SoftVisionXY{AbstractFloat} <: AbstractModel{AbstractFloat}
    T::AbstractFloat
    vision::AbstractFloat
    symmetry::String
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function SoftVisionXY(params)
    @unpack T,vision,symmetry,dt,float_type,rho = params
    T,vision,dt,rho = convert.(float_type,(T,vision,dt,rho))

    if vision ≠ 2π @assert symmetry == "polar" "I am not sure how to interpret a vision cone with nematic symmetry" end

    return SoftVisionXY{float_type}(T,vision,symmetry,dt,zero(float_type),rho)
end

## Propagation Models
mutable struct PropagationForcedXY{AbstractFloat} <: AbstractPropagationModel{AbstractFloat}
    T::AbstractFloat
    Var::AbstractFloat
    symmetry::String
    omegas::AbstractArray{AbstractFloat}
    dt::AbstractFloat
    t::AbstractFloat
    rho::AbstractFloat
end
function PropagationForcedXY(params)
    @unpack L,T,Var,symmetry,dt,float_type,rho  = params
    T,Var,dt,rho = convert.(float_type,(T,Var,dt,rho))

    omegas = zeros(float_type,L) # dummy
    # will be instanciated later during the init, when one knows the type of lattice (1D, 2D)

    rho1 = one(float_type) # for now, we don't want to investigate holes in this model
    return PropagationForcedXY{float_type}(T,Var,symmetry,omegas,dt,zero(float_type),rho1)
end
