## ------------------------ Lattices ------------------------
abstract type AbstractLattice end

mutable struct TriangularLattice <: AbstractLattice
    const L::Int
    periodic::Bool
    const single::Bool
    const metric::String
end
TriangularLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = TriangularLattice(L,periodic,single,metric)

mutable struct SquareLattice <: AbstractLattice
    const L::Int
    periodic::Bool
    const single::Bool
    const metric::String
end
SquareLattice(L::Int;periodic::Bool=true,single::Bool=true,metric::String="euclidian") = SquareLattice(L,periodic,single,metric)

## ------------------------ Functions ------------------------
function dist(lattice::AbstractLattice,pos1,pos2)
    a,b = pos1 ; x,y = pos2
    dx = abs(x-a)
    dy = abs(y-b)
    if lattice.periodic
        dx = min(dx,lattice.L-dx)
        dy = min(dy,lattice.L-dy)
    end
    if     lattice.metric == "euclidian" return sqrt(dx^2 + dy^2)
    elseif lattice.metric == "manhattan" return dx + dy
    elseif lattice.metric == "chebychev" return max(dx,dy)
    else error("Unknown metric !")
    end
end


function distance_matrix(new,old,lattice::AbstractLattice)
    m_new,m_old = length(new),length(old)
    distance_matrix = zeros(m_new,m_old)
    for j in 1:m_old
        for i in 1:m_new
            distance_matrix[i,j] = dist(lattice,new[i],old[j])
        end
    end
    return distance_matrix
end

function offsets(lattice::TriangularLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if even return [(0,1) , (-1,1) , (-1,0) , (0,-1)  , (1,0)  , (1,1)]
    else    return [(0,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1) , (1,0)]
    end
end

function offsets(lattice::SquareLattice,even::Bool)::Vector{Tuple{Int,Int}}
    if lattice.metric in ["manhattan","euclidian"]
        return [(0,1) , (-1,0) , (0,-1) , (1,0)]
    elseif lattice.metric =="chebychev"
        return [(0,1) , (-1,1) , (-1,0) , (-1,-1) , (0,-1) , (1,-1), (1,0) , (1,1)]
    else error("Unknown metric !")
    end
end

function linear_to_square_index(n::Int,L::Int)
    i = div(n-1,L) + 1  # 1 ≤ i ≤ L
    j = mod1(n,L)       # 1 ≤ j ≤ L
    #= formula for reversing the linear indexing. i is the quotient and j
    the reminder of the euclidian division of n by L and the corrections
    deal with the 1-indexing of Julia =#
    return i,j
end

function square_to_linear_index(i::Int,j::Int,L::Int)::Int
    return L*(i-1) + j # 1 ≤ n ≤ L^2
end

function ontheborder(i::Int,j::Int,L::Int)::Bool
    if i == 1 || j == 1 || i == L || j == L
        return true
    else
        return false
    end
end

is_in_bulk(i::Int,j::Int,L::Int)::Bool = !ontheborder(i,j,L)
function at_least_at_distance_X_from_border(i::Int,j::Int,L::Int;X=2)::Bool
    if i < X || j < X || i > L-X+1 || j > L-X+1
        return false
    else
        return true
    end
end

function add_2_positions(pos1::Tuple{T,T},pos2::Tuple{T,T},L::T,should_take_mod::Bool)::Tuple{T,T} where T<:Int
    if should_take_mod return mod1.(pos1 .+ pos2 ,L)
    else return pos1 .+ pos2
    end
end

function mean_2_positions(pos1,pos2,L,should_take_mod::Bool=true)
    a,b = pos1 ; x,y = pos2

    dx = (x - a) #; dx = min(dx,L-dx)
    dy = (y - b) #; dy = min(dy,L-dy)

    if should_take_mod
        if abs(L-dx) < abs(dx) dx = -(L-dx) end
        # if L-dx < dx dx = -(L-dx) end
        if abs(L-dy) < abs(dy) dy = -(L-dy) end
        # if L-dy < dy dy = -(L-dy) end
        estimate_loc_collision = mod1.((a,b) .+ 0.5.*(dx,dy),L)
        # estimate_loc_collision = mod1.(0.5.*(a,b) .+ 0.5.*(x,y),L)
        return estimate_loc_collision
    end
    return (a,b) .+ 0.5.*(dx,dy) # if should_take_mod == false
end
# l = 100
# mean_2_positions((50,50),(60,60),l) == (55,55)
# mean_2_positions((10,10),(90,90),l) == (100,100)
# mean_2_positions((49,66),(51,61),l) == (50.0, 63.5)
