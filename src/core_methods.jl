include("lattices.jl");
include("models.jl");
using StatsBase,Distributions

## ------------------------ Get Neighbours ------------------------
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
        #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
           5  4 (imagine this line is shifhted 1/2 lattice spacing to the left)
        6  X  3
           1  2 (imagine this line is shifhted 1/2 lattice spacing to the left)
        =#
        @inbounds angles =
           # [thetas[i,jp],
           #  thetas[imm,jp],
           #  thetas[imm,j],
           #  thetas[i,jm],
           #  thetas[ip,j],
           #  thetas[ip,jp]]
           [thetas[ip,j],
            thetas[ip,jp],
            thetas[i,jp],
            thetas[imm,jp],
            thetas[imm,j],
            thetas[i,jm]]

    else
        #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
        4  3
        5  X  2 (imagine this line is shifhted 1/2 lattice spacing to the left)
        6  1
        =#
        @inbounds angles =
           # [thetas[i,jp],
           #  thetas[imm,j],
           #  thetas[imm,jm],
           #  thetas[i,jm],
           #  thetas[ip,jm],
           #  thetas[ip,j]]
           [thetas[ip,j],
            thetas[i,jp],
            thetas[imm,j],
            thetas[imm,jm],
            thetas[i,jm],
            thetas[ip,jm]]

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
        #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
           3
        4  X  2
           1
        =#
        @inbounds angles =
            # [thetas[i,jp],
            #  thetas[imm,j],
            #  thetas[i,jm],
            #  thetas[ip,j]]
             [thetas[ip,j],
              thetas[i,jp],
              thetas[imm,j],
              thetas[i,jm]]

    elseif lattice.metric == "chebychev"
        #= Ordering of the nearest neighbours (pure convention, motivated in the file lattices.jl)
        6  5  4
        7  X  3
        8  1  2
        =#
        @inbounds angles =
           # [thetas[i,jp],
           #  thetas[imm,jp],
           #  thetas[imm,j],
           #  thetas[imm,jm],
           #  thetas[i,jm],
           #  thetas[ip,jm],
           #  thetas[ip,j],
           #  thetas[ip,jp]]
            [thetas[ip,j],
            thetas[ip,jp],
            thetas[i,jp],
            thetas[imm,jp],
            thetas[imm,j],
            thetas[imm,jm],
            thetas[i,jm],
            thetas[ip,jm]]
    end
    if model.rho < 1 return filter!(!isnan,angles)
    else return angles
    end
end

function get_neighbours(thetas::Vector{T},model::AbstractPropagationModel{T},lattice::Chain1D,i::Int,bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    L = lattice.L
    if bulk return [thetas[i-1],thetas[i+1]]
    else
        if lattice.periodic return [thetas[mod1(i-1,L)],thetas[mod1(i+1,L)]]
        else
            if i == 1 return [thetas[2]]
            elseif i == L return [thetas[L-1]]
            else return [thetas[i-1],thetas[i+1]]
            end
        end
    end
end

function sum_influence_neighbours(theta::T,i::Int,j::Int,angles_neighbours::Vector{<:T},model::AbstractModel{T},lattice::AbstractLattice)::T where T<:AbstractFloat
    # default case, nothing to worry about
    if isempty(angles_neighbours) return 0.0 # instead of sum(sin,...,init=0) because not possible in Julia 1.3.0 on the cluster I use
    else
        # return sum(sin,sym(model)*(angles_neighbours .- theta)) # 33% times slower
        if model.symmetry == "polar"
            return sum(sin,angles_neighbours .- theta)
        elseif model.symmetry == "nematic"
            return sum(sin,2.0*(angles_neighbours .- theta))
        else error("Symmetry unknown")
        end
    end
end

function sum_influence_neighbours(theta::T,i::Int,j::Int,angles_neighbours::Vector{<:T},model::VisionXY{T},lattice::Abstract2DLattice)::T where T<:AbstractFloat
    weights  = zeros(T,length(angles_neighbours))
    theta0   = mod(theta,modd(model))
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
    if isempty(angles_neighbours) return 0.0 # instead of sum(sin,...,init=0) because not possible in Julia 1.3.0 on the cluster I use
    else return sum(sin.(sym(model)*(angles_neighbours .- theta)) .* weights)
    end
end

function sum_influence_neighbours(theta::T,i::Int,j::Int,angles_neighbours::Vector{<:T},model::SoftVisionXY{T},lattice::Abstract2DLattice)::T where T<:AbstractFloat
    nnn = number_nearest_neighbours(lattice)
    if isempty(angles_neighbours) return 0.0 # instead of sum(sin,...,init=0) because not possible in Julia 1.3.0 on the cluster I use
    else
        ID = ID_projection_angle_onto_lattice(theta,i,j,lattice)
        if model.symmetry == "polar"
            base = sum(sin,angles_neighbours .- theta) * (1. - model.vision)
            correction = nnn*model.vision * sin(angles_neighbours[ID]-theta)
            return base + correction
        # elseif model.symmetry == "nematic"
        #     nnn2 = Int(nnn/2)
        #     base = sum(sin,angles_neighbours .- theta) * (1. - model.vision)
        #     correction_front = nnn2*model.vision * sin(angles_neighbours[ID]-theta)
        #     correction_back  = nnn2*model.vision * sin(angles_neighbours[mod1(ID+nnn2,nnn)]-theta)
        #     return base + correction_front + correction_back
        end
    end
end

## ------------------------ Update Original Models ------------------------
function update!(thetas::AbstractArray{T},model::AbstractModel,lattice::AbstractLattice;tmax) where T<:AbstractFloat
    while model.t < tmax update!(thetas,model,lattice) end
    return thetas
end
update!(thetas::AbstractArray{T},model::AbstractModel,lattice::AbstractLattice,Δt) where T<:AbstractFloat = update!(thetas,model,lattice,tmax=model.t+Δt)

function update!(thetas::Matrix{<:FT},model::Union{LangevinXY{FT},VisionXY{FT},SoftVisionXY{FT}},lattice::Abstract2DLattice) where FT<:AbstractFloat
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
            thetas[i,j] =  θ + dt*sum_influence_neighbours(θ,i,j,angle_neighbours,model,lattice) + sqrt(2T*dt)*randn(FT)
        end
    end

    # On the borders
    ij_in_bulk = false
    if lattice.periodic # if not periodic, juste leave the edge unchanged, otherwise get_neighbours too complicated
        for j in [1,L] , i in 1:L
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*sum_influence_neighbours(θ,i,j,angle_neighbours,model,lattice) + sqrt(2T*dt)*randn(FT)
        end
        for j in 2:L-1 , i in [1,L]
            θ = thetas_old[i,j]
            angle_neighbours = get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            thetas[i,j] =  θ + dt*sum_influence_neighbours(θ,i,j,angle_neighbours,model,lattice) + sqrt(2T*dt)*randn(FT)
        end
    end

    model.t += dt
    return thetas
end

function update!(thetas::Matrix{<:FT},model::MonteCarloXY{FT},lattice::Abstract2DLattice) where FT<:AbstractFloat
    thetas_old = copy(thetas)
    L  = lattice.L
    T  = model.T
    proposals = mod.( 2sqrt(T)*randn(L,L) + thetas ,2pi) # standard deviation two times the standard deviation of the
    # proposals = 2pi*rand(L,L)
    model.symmetry == "polar" ? coeff_symmetry = 1.0 : coeff_symmetry = 2.0
    coeff_symmetry2 = coeff_symmetry / 2.0

    # In Bulk
    ij_in_bulk = true
    for j in 2:L-1
        for i in 2:L-1
            proposal = proposals[i,j] ; theta = thetas[i,j]
            angle_neighbours = 2get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            @fastmath dE = 1.0/coeff_symmetry2 * sin(coeff_symmetry2*(proposal - theta)) * sum(sin,coeff_symmetry2*(proposal + theta .- angle_neighbours))
            if rand() < exp(-dE/T) @inbounds thetas[i,j] = proposal end
        end
    end

    # On the borders
    ij_in_bulk = false
    if lattice.periodic # if not periodic, juste leave the edge unchanged, otherwise get_neighbours too complicated
        for j in [1,L] , i in 1:L
            proposal = proposals[i,j] ; theta = thetas[i,j]
            angle_neighbours = 2get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            @fastmath dE = 1.0/coeff_symmetry2 * sin(coeff_symmetry2*(proposal - theta)) * sum(sin,coeff_symmetry2*(proposal + theta .- angle_neighbours))
            if rand() < exp(-dE/T) @inbounds thetas[i,j] = proposal end
        end
        for j in 2:L-1 , i in [1,L]
            proposal = proposals[i,j] ; theta = thetas[i,j]
            angle_neighbours = 2get_neighbours(thetas_old,model,lattice,i,j,ij_in_bulk)
            @fastmath dE = 1.0/coeff_symmetry2 * sin(coeff_symmetry2*(proposal - theta)) * sum(sin,coeff_symmetry2*(proposal + theta .- angle_neighbours))
            if rand() < exp(-dE/T) @inbounds thetas[i,j] = proposal end
        end
    end

    model.t += 1
    return thetas
end


function update!(thetas::Matrix{<:FT},model::Union{ForcedXY{FT},PropagationForcedXY{FT}},lattice::Abstract2DLattice) where FT<:AbstractFloat
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

function update!(thetas::Matrix{<:FT},model::SPP{FT},lattice::Abstract2DLattice) where FT<:AbstractFloat
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
                collision!(thetas,model,lattice,(i,j),θ,(ic,jc),θc,in_bulk)
            end
        end
    end
    model.t += 1 # here = number of MonteCarlo steps
    return thetas
end

## ------------------------ Update Propagation Models ------------------------
function update!(thetas::Vector{FT},model::AbstractPropagationModel{FT},lattice::Abstract1DLattice)::Vector{FT} where FT<:AbstractFloat
    thetas_old = copy(thetas)
    L  = lattice.L
    dt = model.dt
    T  = model.T

    model.symmetry == "polar" ? symm = 1 : symm = 2

    in_bulk = true
    for i in 2:L-1
        θ = thetas_old[i]
        angle_neighbours = get_neighbours(thetas_old,model,lattice,i,in_bulk)
        thetas[i] =  θ + dt*(model.omegas[i] + sum(sin,symm*(angle_neighbours .- θ)) ) + sqrt(2T*dt)*randn(FT)
    end
    thetas[1] =  thetas_old[1] + dt*(model.omegas[1] + sum(sin,symm*(get_neighbours(thetas_old,model,lattice,1,false) .- thetas_old[1])) ) + sqrt(2T*dt)*randn()
    thetas[L] =  thetas_old[L] + dt*(model.omegas[L] + sum(sin,symm*(get_neighbours(thetas_old,model,lattice,L,false) .- thetas_old[L])) ) + sqrt(2T*dt)*randn()

    model.t += dt
    return thetas
end

## ------------------------ Other Evolution Methods ------------------------

function collision!(thetas::Matrix{<:FT},model::SPP{FT},lattice::AbstractLattice,pos1::Tuple{T,T},theta1::FT,pos2::Tuple{T,T},theta2::FT,bulk::Bool) where {T<:Int,FT<:AbstractFloat}
    # @assert model.symmetry == "nematic" "Energy is only coded for nematic interaction for now !"
    width_proposal = 2sqrt(model.T)
    proposal = width_proposal*randn(FT)+theta1
    i,j = pos1
    if model.algo == "A" # Model A. Align nematically with all NN
        neighbours = 2get_neighbours(thetas,model,lattice,i,j,bulk)
        # dE = -1/2 * ( sum(cos, neighbours .- 2proposal ) - sum(cos, neighbours .- 2theta1 ))
        # dE = sin(proposal - theta1) * sum(sin,proposal + theta1 .- neighbours ) # should be computationally faster
        model.symmetry == "polar" ? coeff_symmetry = 1.0 : coeff_symmetry = 2.0
        coeff_symmetry2 = coeff_symmetry / 2.0
        @fastmath dE = 1.0/coeff_symmetry2 * sin(coeff_symmetry2*(proposal - theta1)) * sum(sin,coeff_symmetry2*(proposal + theta1 .- neighbours))

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

function angle2neighbour(theta::AbstractFloat,i::Int,j::Int,model::AbstractModel,lattice::Abstract2DLattice)
    theta,A = Float64.((theta,model.A)) # Float64.((1,1)) 50x faster than Float64.([1,1])
    direction_motion = direction_of_motion(theta,A)
    NN = project_angle_onto_lattice(direction_motion,i,j,lattice)
    # if model.propulsion == "nematic" && rand(Bool)
    #     NN = (-1) .* NN # choose opposite direction with probability 1/2
    # end
    return add_2_positions((i,j),NN,lattice.L,true) # TODO, check whether false could make it and what runtime gain it would yields
end

function direction_of_motion(theta::T,A::T) where T<:AbstractFloat
    if A == 0 angle = 2π*rand()
    else
        # angle = 1.0/sqrt(A)*randn()+theta  # Wrapped Normal
        # angle = rand(VonMises(theta,A)) # Von Mises, indistinguishable from Wrapped Normal for A > 4
        angle = mod(rand(Cauchy(theta,one(T)/A)),2π) # Wrapped Cauchy, contractile activity

        #= Important Note :
        If instead of centering the variable 'angle' on the variable 'theta',
        one centers it on thetas +π or thetas +π/2, there is no qualitative difference in the movies.
        +1/2 comet-shaped defects are superdiffusive/ballistic, the rest of the behaviour is also left
        unchanged. It basicaly is just a shift in the colours. This means that the details of the coupling
        orientation/polar_propulsion don't matter.
        I take advantage of this fact to plot coherently the field thetas.
        Hence the -π/2 , which is mandatory if one want coherence between the colors
        and the direction of motion once plotted. Because heatmap(thetas') is equivalent to
        a 90° counterclockwise rotation, I choose the neighbour with a 90° clockwise biais (the -pi/2).
        I am aware that this is kind of a weird choice but abandoning the heatmap(thetas') procedure
        complexifies averything else. With this plotting procedure, the visual and the mathematical
        field theta = arctan(y/x) + µ are the same. Otherwise, it seems to me that, for instance,
        µ = 0 is mathematically a source but visually is a counterclockwise vortex etc etc.  =#
    end
    return angle
end
# histogram(mod.(rand(Cauchy(.5,0.25),Int(1E5)),2pi),normalize=true,bins=100)


function project_angle_onto_lattice(angle::AbstractFloat,i::Int,j::Int,lattice::Abstract2DLattice)
    nearest_neighbours = offsets(lattice,iseven(i)) # can be of length 4, 6 or 8
    nb_nn = length(nearest_neighbours)
    index_nearest_neighbour = round(Int,mod(angle,2π)/(2π/nb_nn),RoundNearest) + 1
    if index_nearest_neighbour == nb_nn + 1
        index_nearest_neighbour = 1
    end
    return nearest_neighbours[index_nearest_neighbour]
end

# same as above but only returns the index instead of the offset to reach the said projected neighbour
function ID_projection_angle_onto_lattice(angle::AbstractFloat,i::Int,j::Int,lattice::Abstract2DLattice)::Int
    nearest_neighbours = offsets(lattice,iseven(i)) # can be of length 4, 6 or 8
    nb_nn = length(nearest_neighbours)
    index_nearest_neighbour = round(Int,mod(angle,2π)/(2π/nb_nn),RoundNearest) + 1
    if index_nearest_neighbour == nb_nn + 1
        index_nearest_neighbour = 1
    end
    return index_nearest_neighbour
end

# Meant to relax reconstruction for spotting defects
function relax!(thetas::Matrix{T},model::AbstractModel{T},trelax=0.5) where T<:AbstractFloat
    dummy_dt = T(1E-2)
    dummy_model = LangevinXY{T}(zero(T),model.symmetry,dummy_dt,zero(T),model.rho)
    dummy_lattice = SquareLattice(size(thetas,1),true,true,"chebychev")
    update!(thetas,dummy_model,dummy_lattice,trelax)
end
