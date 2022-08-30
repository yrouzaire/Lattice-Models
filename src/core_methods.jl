include("lattices.jl");
include("models.jl");

## ------------------------ Get Neighbours ------------------------
# NOTE : no need to define get_neighbours(model::AbstractModel,lattice::AbstractLattice ...)
function get_neighbours(thetas::Matrix{<:T},model::AbstractModel{T},lattice::TriangularLattice,i::Int,j::Int,bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    L = lattice.L

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

function get_neighbours(thetas::Matrix{<:T},model::AbstractModel{T},lattice::SquareLattice,i::Int,j::Int,bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    L = lattice.L

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
            thetas[imm,j],
            thetas[i,jm],
            thetas[ip,j]]

    else
        @inbounds angles =
           [thetas[i,jp],
            thetas[imm,j],
            thetas[i,jm],
            thetas[ip,j]]
    end
    return angles
end

function sum_influence_neighbours(theta::T,angles_neighbours::Vector{<:T},model::AbstractModel{T},lattice::AbstractLattice)::T where T<:AbstractFloat
    # default case, nothing to worry about
    return sum(sin,angles_neighbours .- theta,init=0) # init = 0 in case angles_neighbours is empty
end

function sum_influence_neighbours(theta::T,angles_neighbours::Vector{<:T},model::VisionXY{T},lattice::AbstractLattice)::T where T<:AbstractFloat
    weights  = zeros(T,length(angles_neighbours))
    symm     = sym(model)
    theta0   = mod(theta,symm)
    if     isa(lattice,TriangularLattice) constant = T(π/3)
    elseif isa(lattice,SquareLattice)     constant = T(π/2)
    end

    for n in 1:length(angles_neighbours)
        dtheta     = theta0 - (n-1)*constant
        dtheta_abs = abs(dtheta)
        arclengt   = min(symm-dtheta_abs,dtheta_abs)

        if arclengt ≤ model.vision/2
            @inbounds weights[n] = 1.0
        end
    end

    return sum(sin.(angles_neighbours .- theta) .* weights , init = 0) # init = 0 in case angles_neighbours is empty
end


## ------------------------ Update ------------------------
# NOTE : no need to define update!(model::AbstractModel,lattice::AbstractLattice)
function update!(thetas::Matrix{<:T},model::AbstractModel,lattice::AbstractLattice,tmax::Number) where T<:AbstractFloat
    while model.t < tmax
        update!(thetas,model,lattice)
    end
    return thetas
end

function update!(thetas::Matrix{<:FT},model::Union{XY{FT},VisionXY{FT}},lattice::AbstractLattice) where FT<:AbstractFloat
    thetas_old = copy(thetas)
    L  = lattice.L
    dt = model.dt
    T  = model.T
    # In Bulk
    ij_in_bulk = true
    for j in 2:L-1
        for i in 2:L-1
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*sum_influence_neighbours(θ,angle_neighbours,model,lattice) + sqrt(2T*dt)*randn(FT)
        end
    end

    # On the borders
    ij_in_bulk = false
    if lattice.periodic # if not periodic, juste leave the edge unchanged, otherwise get_neighbours too complicated
        for j in [1,L] , i in 1:L
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*sum_influence_neighbours(θ,angle_neighbours,model,lattice) + sqrt(2T*dt)*randn(FT)
        end
        for j in 2:L-1 , i in [1,L]
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*sum_influence_neighbours(θ,angle_neighbours,model,lattice) + sqrt(2T*dt)*randn(FT)
        end
    end

    model.t += dt
    return thetas
end

function update!(thetas::Matrix{<:FT},model::AXY{FT},lattice::AbstractLattice) where FT<:AbstractFloat
    thetas_old = copy(thetas)
    L  = lattice.L
    dt = model.dt
    T  = model.T
    # In Bulk
    ij_in_bulk = true
    for j in 2:L-1
        for i in 2:L-1
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*(model.omegas[i,j] + sum(sin,angle_neighbours .- θ)) + sqrt(2T*dt)*randn(FT)
        end
    end

    # On the borders
    ij_in_bulk = false
    if lattice.periodic # if not periodic, juste leave the edge unchanged, otherwise get_neighbours too complicated
        for j in [1,L] , i in 1:L
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*(model.omegas[i,j] + sum(sin,angle_neighbours .- θ)) + sqrt(2T*dt)*randn(FT)
        end
        for j in 2:L-1 , i in [1,L]
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*(model.omegas[i,j] + sum(sin,angle_neighbours .- θ)) + sqrt(2T*dt)*randn(FT)
        end
    end

    model.t += dt
    return thetas
end
