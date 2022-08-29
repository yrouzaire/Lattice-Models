include("lattices.jl");
include("models.jl");

## ------------------------ Get Neighbours ------------------------
# NOTE : no need to define get_neighbours(model::AbstractModel,lattice::AbstractLattice ...)
# TODO deal with non periodic BC
function get_neighbours(model::AbstractModel{T},lattice::TriangularLattice,i::Int,j::Int,bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    L = model.L
    thetas = model.thetas

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

function get_neighbours(model::VisionXY{T},lattice::TriangularLattice,i::Int,j::Int,bulk::Bool=false)::Vector{T} where T<:AbstractFloat
    all_angles = invoke(get_neighbours, Tuple{AbstractModel{T},typeof(lattice),Int,Int,Bool}, model,lattice,i,j,bulk)
    neighbours_in_vision_cone = T[]
    symm = sym(model)
    theta0 = mod(model.thetas[i,j],symm)
    for n in 1:length(all_angles)
        dtheta = theta0 - (n-1)*π/3 # because TriangularLattice
        dtheta_abs = abs(dtheta)
        arcleng = min(symm-dtheta_abs,dtheta_abs)

        if arcleng ≤ model.vision/2
            push!(neighbours_in_vision_cone,all_angles[n])
        end
    end
    return neighbours_in_vision_cone
    #= Small note on the invoke function used above. To avoid infinite loops
    of get_neighbours(model::VisionXY{T},,lattice::AbstractLattice ....) calling
    itself over and over, I forced it to invoke the more general one, namely
    get_neighbours(model::AbstratcModel{T},lattice::AbstractLattice ....) .
    I then refine the general result to satisfy the vision cone.

    Hereafter, a small working example of the use of the invoke function.
    ``
    ft(x::Int) = "Int with type(x) = $(typeof(x))"
    ft(x:: Any) = "Any with type(x) = $(typeof(x))"
    a = 1 ; b = "yy"
    ft(a)
    ft(b)
    invoke(ft,Tuple{Any},a)
    ``
    =#

    #= Another Note : I tried to be more general by defining
    get_neighbours(model::VisionXY{T},lattice::AbstractLattice) but error
    while running :
    MethodError: get_neighbours(::VisionXY{Float32}, ::TriangularLattice ... ) is ambiguous. Candidates:
    get_neighbours(model::AbstractModel{T}, lattice::TriangularLattice ...)
    get_neighbours(model::VisionXY{T}, lattice::AbstractLattice, ...)

    Possible fix, define get_neighbours(::VisionXY{T}, ::TriangularLattice, ...)
    =#
end


## ------------------------ Update ------------------------
# NOTE : no need to define update!(model::AbstractModel,lattice::AbstractLattice)
function update!(model::Union{AXY{FT},XY{FT}},lattice::AbstractLattice) where FT<:AbstractFloat
    # In Bulk
    ij_in_bulk = true
    for j in 2:model.L-1
        for i in 2:model.L-1
            langevin_update!(model,i,j,ij_in_bulk)
        end
    end
    # On the borders
    ij_in_bulk = false
    for j in [1,L] , i in 1:model.L
        langevin_update!(model,i,j,ij_in_bulk)
    end
    for j in 2:model.L-1 , i in [1,L]
        langevin_update!(model,i,j,ij_in_bulk)
    end

    model.thetas = model.thetas_new
    return thetas
end

function langevin_update!(model::XY{FT},i::Int,j::Int,ij_in_bulk::Bool) where FT<:AbstractFloat
    θ = model.thetas[i,j]
    angle_neighbours = get_neighbours(model,i,j,ij_in_bulk)
    model.thetas_new[i,j] =  θ + model.dt*sum(sin,angle_neighbours .- θ) + sqrt(2*model.T*model.dt)*randn(FT)
    return nothing
end

function langevin_update!(model::AXY{FT},i::Int,j::Int,ij_in_bulk::Bool) where FT<:AbstractFloat
    θ = model.thetas[i,j]
    angle_neighbours = get_neighbours(model,i,j,ij_in_bulk)
    model.thetas_new[i,j] =  θ + model.dt*(model.omegas[i,j]  + sum(sin,angle_neighbours .- θ)) + sqrt(2*model.T*model.dt)*randn(FT)
    return nothing
end
