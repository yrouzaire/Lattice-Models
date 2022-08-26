include("lattice_general.jl");
include("models.jl");

## ------------------------ Generic methods ------------------------
function get_neighbours_triangular(thetas::Matrix{T1},L::T2,i::T2,j::T2,bulk::Bool)::Vector{<:T1} where {T1<:AbstractFloat,T2<:Int}
    # convention depuis la droite et sens trigo
    if bulk
        jm  = j-1
        jp  = j+1
        imm = i-1
        ip  = i+1
    else
        jm  = mod1(j-1,L)
        jp  = mod1(j+1,L)
        imm = mod1(i-1,L)
        ip  = mod1(i+1,L)
    end

    if iseven(i)
        @inbounds angles =
           [thetas[i,jp],
            thetas[imm,jp],
            thetas[imm,j],
            thetas[i,jm],
            thetas[ip,j],
            thetas[ip,jp]]

    else
        @inbounds angles =
           [thetas[i,jp],
            thetas[imm,j],
            thetas[imm,jm],
            thetas[i,jm],
            thetas[ip,jm],
            thetas[ip,j]]
    end
    return angles
end

function get_neighbours(model::AbstractModel,lattice::TriangularLattice,i::Int,j::Int;bulk::Bool=false)::Vector{<:AbstractFloat}
    # The method called by default for all TriangularLattices if the model is not VisionXY
    @assert lattice.periodic  #TODO
    return get_neighbours_triangular(model.thetas,model.L,i,j,bulk)
end

function update!(model::AbstractModel,lattice::TriangularLattice)

end

## ---------------------- Specific methods for VisionXY ----------------------
function get_neighbours(model::VisionXY{T},lattice::TriangularLattice,i::Int,j::Int;bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    # The method called if and only if the model is VisionXY
    @assert lattice.periodic  #TODO
    angles = get_neighbours_triangular(model.thetas,model.L,i,j,bulk)
    theta = mod(model.thetas[i,j],sym(model))
    result = T[]
    for i in 1:length(angles)
        dtheta = theta - (i-1)*π/3
        dtheta_abs = abs(dtheta)
        arcleng = min(sym(model)-dtheta_abs,dtheta_abs)

        if arcleng ≤ model.vision/2
            push!(result,angles[i])
        end
    end
    return result
end
