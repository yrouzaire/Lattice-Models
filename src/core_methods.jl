include("lattices.jl");
include("models.jl");
using StatsBase,Distributions

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
    if model.rho < 1 return filter!(!isnan,angles)
    else return angles
    end
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

    if lattice.metric in ["manhattan","euclidian"]
        @inbounds angles =
           [thetas[i,jp],
            thetas[imm,j],
            thetas[i,jm],
            thetas[ip,j]]

    elseif lattice.metric == "chebychev"
        @inbounds angles =
           [thetas[i,jp],
            thetas[imm,jp],
            thetas[imm,j],
            thetas[imm,jm],
            thetas[i,jm],
            thetas[ip,jm],
            thetas[ip,j],
            thetas[ip,jp]]
    end
    if model.rho < 1 return filter!(!isnan,angles)
    else return angles
    end
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
    while model.t < tmax update!(thetas,model,lattice) end
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

# Meant to relax reconstruction for spotting defects
function relax!(thetas::Matrix{T},model::AbstractModel{T}) where T<:AbstractFloat
    dummy_dt = T(1E-2)
    trelax = T(1.0)
    dummy_model = XY{T}(zero(T),model.symmetry,dummy_dt,zero(T),model.rho)
    dummy_lattice = SquareLattice(size(thetas,1),true,true,"euclidian")
    update!(thetas,dummy_model,dummy_lattice,trelax)
end

function update!(thetas::Matrix{<:FT},model::ForcedXY{FT},lattice::AbstractLattice) where FT<:AbstractFloat
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

function update!(thetas::Matrix{<:FT},model::MovingXY{FT},lattice::AbstractLattice) where FT<:AbstractFloat
    L = lattice.L
    order_trials = StatsBase.shuffle(1:L^2)
    for n in order_trials
        i,j = linear_to_square_index(n,L)
        in_bulk = 2 < i < L-1 &&  2 < j < L-1 # si "seulement" 1 < x,y < L , bug pour avoir le voisin d'un voisin

        θ = thetas[i,j]
        if !isnan(θ) && ( in_bulk || lattice.periodic) # translation : if on the border AND lattice not periodic, do nothing

            ic,jc = angle2neighbour(θ,i,j,model,lattice)
            θc = thetas[ic,jc]

            if isnan(θc) # destination empty, let's go
                thetas[ic,jc] = θ
                thetas[i,j] = NaN
            else # destination occupied
                collision!(thetas,model,(i,j),θ,(ic,jc),θc,in_bulk)
            end
        end
    end
    model.t += 1 # here = number of MonteCarlo steps
    return thetas
end

function collision!(thetas::Matrix{<:FT},model::MovingXY{FT},pos1::Tuple{T,T},theta1::FT,pos2::Tuple{T,T},theta2::FT,bulk::Bool) where {T<:Int,FT<:AbstractFloat}
    @assert model.symmetry == "nematic" "Energy is only coded for nematic interaction for now !"
    proposal = model.width_proposal*randn(FT)+theta1
    i,j = pos1
    if model.algo == "A" # Model A. Align nematically with all NN
        neighbours = 2get_neighbours(thetas,model,lattice,i,j,bulk)
        dE = -1/2 * ( sum(cos, neighbours .- 2proposal ) - sum(cos, neighbours .- 2theta1 ))
        # dE = sin(proposal - theta1) * sum(sin,proposal + theta1 .- neighbours ) # should be computationally faster
    elseif model.algo == "B" # align nematically wrt collided spin
        dE = -1/2 * ( cos(2(theta2 - proposal)) - cos(2(theta2 - theta1)))
    elseif model.algo == "C" # align F/AF wrt collided spin
        ccos = cos(theta1 - theta2)
        if ccos > 0  J = +1.0 # ferro
        else J = -1.0 # antiferro
        end
        dE = -J*(cos(theta1 - proposal) - ccos)
    elseif model.algo == "CA"
        ccos = cos(theta1 - theta2)
        if ccos > 0  J = +1.0 # ferro
        else J = -1.0 # antiferro
        end
        dE = -J*(cos(theta1 - proposal) - ccos)
        if rand() < exp(-dE/model.T) # Metropolis Monte Carlo
            @inbounds thetas[i,j] = proposal
        end
        proposal = model.width_proposal*randn(FT)+thetas[i,j] # new proposal
        neighbours = 2get_neighbours(thetas,model,lattice,i,j,bulk)
        dE = -1/2 * ( sum(cos, neighbours .- 2proposal ) - sum(cos, neighbours .- 2theta1 ))
        # the last Monte Carlo step, corresponding to the case "A", is performed by the last step (common to all cases)
    else error("Unknown Model")
    end
    if rand() < exp(-dE/model.T) # Metropolis Monte Carlo
        @inbounds thetas[i,j] = proposal
    end
    return thetas
end

function angle2neighbour(theta::AbstractFloat,i::Int,j::Int,model::AbstractModel,lattice::AbstractLattice)
    theta,A = Float64.((theta,model.A)) # Float64.((1,1)) 50x faster than Float64.([1,1])
    direction_motion = direction_of_motion(theta,A)
    NN = project_angle_onto_lattice(direction_motion,i,j,lattice)
    if model.propulsion == "nematic" && rand(Bool)
        NN = (-1) .* NN # choose opposite direction with probability 1/2
    end
    return add_2_positions((i,j),NN,lattice.L,true) # TODO, check whether false could make it and what runtime gain it would yields
end

function direction_of_motion(theta::T,A::T) where T<:AbstractFloat
    if A == 0 angle = 2π*rand()
    else
        # angle = 1.0/sqrt(A)*randn()+theta  # Wrapped Normal
        angle = rand(VonMises(theta,A)) # Von Mises, indistinguishable from Wrapped Normal for A > 4
    end
    return angle
end


function project_angle_onto_lattice(angle::AbstractFloat,i::Int,j::Int,lattice::AbstractLattice)
    nearest_neighbours = offsets(lattice,iseven(i)) # can be of length 4, 6 or 8
    nb_nn = length(nearest_neighbours)
    index_nearest_neighbour = round(Int,mod(angle,2π)/(2π/nb_nn),RoundNearest) + 1
    if index_nearest_neighbour == nb_nn + 1
        index_nearest_neighbour = 1
    end
    return nearest_neighbours[index_nearest_neighbour]
end
